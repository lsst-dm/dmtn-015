# Copyright (c) 2012-2016 by the GalSim developers team on GitHub
# https://github.com/GalSim-developers
#
# This file a modified copy of part of
# GalSim: The modular galaxy image simulation toolkit.
# https://github.com/GalSim-developers/GalSim
#
# GalSim is free software: redistribution and use in source and binary forms,
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions, and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions, and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.
#


import numpy
import scipy.stats
import matplotlib.pyplot
import galsim

numpy.random.seed(500)

N_SIGMA_DEPTH = 5    # use 5-sigma (point source) as the measure of depth

# Kwargs to pass to drawImage when creating PSF images.
PSF_DRAW_KWARGS = dict(nx=20, ny=20, scale=0.2, method='no_pixel')


class CoaddMocker(object):
    """Interface specification for classes that know how to mock up an effective
    PSF and variance for a particular coaddition algorithm.
    """

    def mockCoadd(self, variances, fwhms, psfs):
        """Compute the effective PSF and per-pixel variance of a coadd.

        This method must be implemented by subclasses to define the coadd
        algorithm.
        """
        raise NotImplementedError()

    def selectInputs(self, variances, fwhms):
        """Return a mask object (boolean array or slice) that selects images
        that should go into a coadd, given the FWHMs and variances of the input
        exposures.

        Default implementation returns a slice that selects all images.
        """
        slice(None)


class DirectCoaddMocker(object):
    """Coadd mocker for direct coaddition: simply adding images with some weight.
    """

    def __init__(self, included_fraction=1.0):
        self.included_fraction = included_fraction

    def mockCoadd(self, variances, fwhms, psfs):
        weights = 1.0 / (fwhms**2 * variances)
        weights /= weights.sum()
        coadd_psf = galsim.Sum([psf*weight for psf, weight in zip(psfs, weights)])
        coadd_variance = (variances*weights*weights).sum()
        return coadd_variance, coadd_psf

    def selectInputs(self, variances, fwhms):
        cutoff = numpy.percentile(fwhms, 100*self.included_fraction)
        return fwhms < cutoff


class PSFMatchedCoaddMocker(object):
    """Coadd mocker for PSF-matched coaddition, in which each exposure's
    PSF is convolved with a kernel that matches it to the worst input PSF.
    """


    def __init__(self, included_fraction=1.0):
        self.included_fraction = included_fraction

    def mockCoadd(self, variances, fwhms, psfs):
        weights = 1.0 / (fwhms**2 * variances)
        weights /= weights.sum()
        coadd_psf = psfs[fwhms.argmax()]
        # We ignore transfer from variance to covariance from convolution
        # with the matching kernel, because this toy model of coaddition
        # doesn't track covariance at all.  Treating it all as variance
        # should be a better approximation than calculating the transfer
        # and then throwing away the variance as long as the matching
        # kernels are typically smaller than the PSFs.
        coadd_variance = (variances*weights*weights).sum()
        return coadd_variance, coadd_psf

    def selectInputs(self, variances, fwhms):
        cutoff = numpy.percentile(fwhms, 100*self.included_fraction)
        return fwhms < cutoff


class KaiserCoaddMocker(object):
    """Coadd mocker for optimal coaddition in Fourier space;
    see http://adsabs.harvard.edu/abs/2015arXiv151206879Z.
    """

    def mockCoadd(self, variances, fwhms, psfs):
        weights = 1.0 / variances
        weights /= weights.sum()
        coadd_variance = 1.0 / (1.0 / variances).sum()
        coadd_psf = galsim.FourierSqrt(
            galsim.Sum([
                galsim.AutoCorrelate(psf)*weight
                for psf, weight in zip(psfs, weights)
            ])
        )
        return coadd_variance, coadd_psf

    def selectInputs(self, variances, fwhms):
        return slice(None, None, None)


class CoaddMetricCalculator(object):
    """Object that runs several coadd mockers on a set of input PSFs and variances
    representing a single point on the sky, computing the effective FWHM of the
    coadd PSF (from its effective area) and the variance in the coadd pixels.
    """

    def __init__(self, mockers=None, included_fraction=1.0):
        if mockers is None:
            mockers = {
                "direct": DirectCoaddMocker(included_fraction),
                "psf-matched": PSFMatchedCoaddMocker(included_fraction),
                "kaiser": KaiserCoaddMocker(),
            }
        self.mockers = dict(mockers)

    def makePSF(self, fwhm):
        """Create a galsim.GSObject that represents a PSF with the given FWHM.

        For simplicity and speed we just use Gaussian.
        """
        return galsim.Gaussian(fwhm=fwhm)

    def buildInputs(self, fwhms, depths):
        """Build input PSFs and variances for a set of input images from their
        FWHMs and 5-sigma magnitude limits.
        """
        psfs = numpy.zeros(len(fwhms), dtype=object)
        n_effs = numpy.zeros(len(fwhms), dtype=float)
        for n, (fwhm, depth) in enumerate(zip(fwhms, depths)):
            psfs[n] = self.makePSF(fwhm)
            psf_image = psfs[n].drawImage(**PSF_DRAW_KWARGS)
            n_effs[n] = psf_image.array.sum()**2 / (psf_image.array**2).sum()
        flux_limit = 10**(-0.4*depths)
        variances = (flux_limit/N_SIGMA_DEPTH)**2 / n_effs
        fwhm_factor = numpy.median(fwhms / n_effs**0.5)
        return psfs, variances, fwhm_factor

    def computeMetrics(self, coadd_psf, coadd_variance, fwhm_factor):
        """Given a coadd PSF (GSObject), per-pixel variance, compute the
        effective FWHM and 5-sigma magnitude limit of the coadd.

        Effective FWHM is computed as a simple scaling of the square root of
        the PSF effective area; effective area is a more meaningful measure of
        PSF size, but FWHM is more readily understood by humans.  The scaling
        factor for this conversion is given by the fwhm_factor argument.
        """
        coadd_psf_image = coadd_psf.drawImage(**PSF_DRAW_KWARGS)
        coadd_n_eff = coadd_psf_image.array.sum()**2 / (coadd_psf_image.array**2).sum()
        coadd_depth = -2.5*numpy.log10(
            N_SIGMA_DEPTH * (coadd_variance*coadd_n_eff)**0.5
        )
        return fwhm_factor*coadd_n_eff**0.5, coadd_depth

    def __call__(self, fwhms, depths):
        """Compute coadd PSF metrics for input exposures defined by the given
        PSF FWHMs and 5-sigma magnitude limits.

        Returns a dict of metrics with keys "<name>.fwhm" and "<name>.depth",
        along with the fwhm_factor used to compute effective FWHM from PSF
        effective area.
        """
        psfs, variances, fwhm_factor = self.buildInputs(fwhms, depths)
        result = {}
        for name, mocker in self.mockers.iteritems():
            mask = mocker.selectInputs(variances, fwhms)
            coadd_variance, coadd_psf = mocker.mockCoadd(variances[mask], fwhms[mask], psfs[mask])
            coadd_fwhm, coadd_depth = self.computeMetrics(coadd_psf, coadd_variance,
                                                          fwhm_factor)
            result["{}.fwhm".format(name)] = coadd_fwhm
            result["{}.depth".format(name)] = coadd_depth
        return result, fwhm_factor


def compareCoadds(depth=24.7, depth_scatter=0.2, fwhm=0.7, fwhm_scatter=0.2,
                  n_exposures=200, n_realizations=100):
    """Generate several realizations of mock inputs and plot histograms of coadd quality.
    """
    depths = scipy.stats.norm(depth, depth_scatter).rvs(size=(n_realizations, n_exposures))
    fwhms = scipy.stats.lognorm(s=fwhm_scatter, scale=fwhm).rvs(size=(n_realizations, n_exposures))

    result = {}

    calc = CoaddMetricCalculator(
        mockers={
            "PSF-matched, best 80%": PSFMatchedCoaddMocker(0.8),
            "PSF-matched, all": PSFMatchedCoaddMocker(1.0),
            "Direct, best 80%": DirectCoaddMocker(0.8),
            "Direct, all": DirectCoaddMocker(1.0),
            "Kaiser": KaiserCoaddMocker(),
        }
    )
    for name in calc.mockers:
        result["{}.fwhm".format(name)] = numpy.zeros(n_realizations, dtype=float)
        result["{}.depth".format(name)] = numpy.zeros(n_realizations, dtype=float)
    for n in xrange(n_realizations):
        local, fwhm_factor = calc(fwhms[n], depths[n])
        for k, v in local.iteritems():
            result[k][n] = local[k]

    plot_kwargs = [
        ("PSF-matched, best 80%", dict(histtype='step', edgecolor='b')),
        ("PSF-matched, all", dict(histtype='stepfilled', facecolor='b', linewidth=0, alpha=0.5)),
        ("Direct, best 80%", dict(histtype='step', edgecolor='r')),
        ("Direct, all", dict(histtype='stepfilled', facecolor='r', linewidth=0, alpha=0.5)),
        ("Kaiser", dict(histtype='stepfilled', facecolor='g', linewidth=0, alpha=0.5)),
    ]

    input_kwargs = dict(histtype='stepfilled', facecolor='c', linewidth=0, alpha=0.5, bins=150, normed=True)

    for k, d in plot_kwargs:
        d.update(bins=150, normed=True)

    fig = matplotlib.pyplot.figure(figsize=(12,12))
    ax1 = fig.add_subplot(2, 1, 1)
    for name, kwargs in plot_kwargs:
        ax1.hist(result["{}.depth".format(name)], label=name, range=(24.0, 28.0), **kwargs)
    ax1.hist(depths.ravel(), label="Input Exposures", range=(24.0, 28.0), **input_kwargs)
    ax1.set_xlabel("5-sigma magnitude limit")
    ax1.set_xlim(24.0, 28.0)
    ax1.set_ylim(0, 25)

    ax2 = fig.add_subplot(2, 1, 2)
    for name, kwargs in plot_kwargs:
        ax2.hist(result["{}.fwhm".format(name)], label=name, range=(0.5, 1.2), **kwargs)
    ax2.hist(fwhms.ravel(), label="Input Exposures", range=(0.5, 1.2), **input_kwargs)
    ax2.set_xlabel("PSF Effective FWHM")
    ax2.set_xlim(0.5, 1.2)
    ax2.set_ylim(0, 50)
    ax2.legend()

    fig.savefig("_static/comparison.png")
    matplotlib.pyplot.show()

    return result, fig

if __name__ == "__main__":
    compareCoadds()
:tocdepth: 2


Notation
========

Given multiple overlapping images :math:`\{z_1({\bf r}), z_2({\bf r}), ...\}`, we designate by :math:`R({\bf r})` the set of image indices that overlap a point on the sky :math:`{\bf r}`.  Each image is related to the true sky :math:`f({\bf r})` via its point spread function (PSF) :math:`\phi_i({\bf r}, {\bf s})` by

.. math::
  z_i({\bf r}) = \int \! \phi_i({\bf r}, {\bf s}) \, f({\bf s}) \, d{\bf s}

We assume that all images have already been resampled to a common coordinate system, and have Gaussian noise described by a covariance matrix :math:`C_i({\bf r}, {\bf s})`.  Because astronomical images typically have approximately uncorrelated noise only in their original coordinate system, we should not assume that these covariance matrices are diagonal.


Lossy Coaddition Algorithms
===========================

Direct Coadds
-------------

The simplest way to build a coadd is to simply form a linear combination of the images that overlap a point on the sky.  The coadd image :math:`z_\mathrm{dir}({\bf r})` is then

.. math::
  z_\mathrm{dir}({\bf r}) = \sum_{i \in R({\bf r})} w_i({\bf r}) \, z_i({\bf r})

where :math:`w_i` is an arbitrary weight for each image, defined to sum to one at each point.  We cannot assume the weights are constant over a single input image; this would make it impossible for the weights to sum to one on both sides of a boundary in which the number of input images changes.  It is equally important that the weight function vary only on scales much larger than the size of the PSF; without this it the coadd has no well-defined effective PSF.  The weight can be chosen to optimize some measure of signal-to-noise ratio (SNR) on the coadd, popular choices include

.. math::
  w_i \propto \frac{1}{\sigma_i^2}

to optimize the average per-pixel SNR and

.. math::
  w_i \propto \frac{|\phi_i|^2}{\sigma_i^2}

to optimize the SNR of a point source.  In both cases, :math:`\sigma^2` is some spatially-averaged scalar reduction of :math:`C_i`, while :math:`|\phi_i|^2` is the inner product of the average PSF (and the inverse of its effective area).

If the weights are proportional to exposure time and the input images are observed back-to-back, the direct coadd is mathematically equivalent to a single longer observation in the limit of perfectly linear detectors.

The effective PSF on the coadd and pixel covariance matrix are simple to compute:

.. math::
  \phi_\mathrm{dir}({\bf r},{\bf s}) =
    \sum_{i \in R({\bf r})} w_i({\bf r}) \, \phi_i({\bf r}, {\bf s})
  :label: eq:coaddpsf

.. math::
  C_\mathrm{dir}({\bf r}, {\bf s}) =
    \sum_{i \in \left[R({\bf r}) \cap R({\bf s})\right]} \!\!\!
        w_i({\bf r}) \, w_i({\bf s})
        \, C_i({\bf r}, {\bf s})

These quantities would be many times larger than the coadd itself if evaluated on every pixel, making direct evaluation impractical.  They are discontinuous at the boundaries of input images (and masked regions within them), making interpolation problematic as well.

The solution we have adopted for PSF models has been referred as both :ref:`CoaddPsf <coaddpsf>` (from :py:class:`lsst.meas.algorithms.CoaddPsf`) and :ref:`StackFit <stackfit>` (after the shear estimation technique where it was developed), and is essentially a form of lazy evaluation.  When a PSF model image is requested at a point, we simply evaluate the PSF models for all of the input images at that point, transform them to the correct coordinate system, and compute the weighted sum on demand.  We typically assume that the PSF is constant over the scale of a single astronomical object, and hence this reduces the number of PSF model evaluations from the number of pixels to the number of detected objects.  When an object lies on a boundary or a region with masked pixels, the true PSF is discontinuous and the constant-PSF assumption is not valid.  At present, we simply flag objects for which this is true, but this may not work when the number of input images is large; in this regime the number of border and masked regions increases, though the severity of the discontinuities decreases as well.

Our current approach for coadd uncertainty propagation is to compute and store only the variance.  We will likely expand this in the future to storing some approximation to the covariance (e.g. by modeling it as constant within regions where the number of input images is constant).

Direct coadds are lossy, requiring some trade-off between image quality (PSF size) and depth (SNR).  This can be easily seen from :eq:`eq:coaddpsf`: including an image with a PSF larger than the current weighted mean PSF always increases the size of the final PSF, regardless of the depth of the new image.


PSF-Matched Coadds
------------------


Outlier Rejection and Nonlinear Statistics
------------------------------------------


Exact Coaddition Algorithms
===========================

Likelihood Coadds
-----------------

Decorrelated Coadds
-------------------

Kaiser Coadds
-------------

Constant PSF Coadds
-------------------


Coadds for Source Detection
===========================

Detection Maps
--------------

Chi-Squared Coadds
--------------------



Glossary
========

.. _chisq_coadd:

Chi-Squared Coadd
  A cross-band coadd that is designed for detecting objects by rejecting the null hypothesis that a pixel contains only sky.  See [Szalay1999]_.

.. _coaddpsf:

CoaddPsf
  A procedure for generating the PSF model at a point on a direct coadd by lazily evaluating the PSF models of the input at that point, then warping and combining them with the same weights used to build the coadd itself.  Originally developed by [Jee2011]_ as part of :ref:`StackFit <stackfit>`.

.. _constant_psf_coadd:

Constant-PSF Coadd
  Any coadd that has been designed to have a constant (spatially non-variable).  This includes :ref:`PSF-matched coadds <psf_matched_coadd>`, but we will frequently use this term instead as shorthand for a partially :ref:`decorrelated coadd <decorrelated_coadd>` with a constant PSF, in which the noise in a :ref:`likelihood coadd <likelihood_coadd>` is only partially decorrelated in order to produce an image with a constant PSF.  A :ref:`Kaiser coadd <kaiser_coadd>` is technically such a coadd, but only because it assumes constant input PSFs.

.. _decorrelated_coadd:

Decorrelated Coadd
  An optimal coadd produced by decorrelating the noise in a :ref:`Likelihood Coadd <likelihood_coadd>`.  The :ref:`Kaiser Coadd <kaiser_coadd>` is a special case that relies restrictive assumptions about the input; the general algorithm can be described mathematically but is computationally impractical without some other approximation.

.. _deep_coadd:

Deep Coadd
  A lossy coadd produced using all but the very worst-seeing images.  Contrast with :ref:`Good-Seeing Coadd <good_seeing_coadd>`.

.. _detection_map:

Detection Map
  An image that can be thresholded to detect sources under the assumption that they are unblended point sources, formed by convolving an image by the transpose of its PSF and dividing each pixel by its variance.  It can also be built by dividing a :ref:`likelihood coadd <likelihood_coadd>` by its variance.

.. _direct_coadd:

Direct Coadd
  A lossy coadd built as a linear combination of images with no change to their PSFs.  If the weights are just the exposure times of the image, this is (locally) equivalent to a single long exposure.  The PSF of a direct coadd is discontinuous at the boundaries of input images, requiring an approach like :ref:`CoaddPsf <coaddpsf>` to model it.  This coadd is lossy, requiring some tradeoff to be made (in selecting inputimages) between depth and image quality.  Noise in a direct coadd is correlated only by image resampling.

.. _good_seeing_coadd:

Good-Seeing Coadd
  A lossy coadd produced using only input images with good seeing.  Constrast with "Deep Coadd."

.. _kaiser_coadd:

Kaiser Coadd
  An optimal coadd built by decorrelating a :ref:`Likelihood Coadd <likelihood_coadd>` after assuming input images have uncorrelated white noise, constant PSFs, and no missing pixels or boundaries.  Origin is [Kaiser2001]_, an unpublished Pan-STARRS white paper.  Special case of :ref:`Decorrelated Coadd <decorrelated_coadd>`.

.. _likelihood_coadd:

Likelihood Coadd
  An optimal coadd built as a linear combination of images that have been convolved with the transpose of their PSFs.  This procedure correlates noise, but the resulting image is optimal for isolated point source detection even if only the variance is propagated and stored (see :ref:`Detection Map <detection_map>`).  For other applications (including producing :ref:`Decorrelated Coadds <decorrelated_coadd>`), the full covariance must be propagated.

.. _multifit:

MultiFit
  An approach to source measurement (especially weak lensing shear estimation) that fits the same model to all input images directly, after transforming the model to the coordinate system of each image and convolving with that image's PSF.  Formally optimal (for valid models).  Contrast with :ref:`StackFit <stackfit>`.

.. _proper_image:

Proper Image
  An image with uncorrelated white noise; see [Zackay2015]_.

.. _psf_matched_coadd:

PSF-Matched Coadd
  A lossy coadd built by combining images only after they have been reconvolved to a common, constant PSF.  This either degrades all images to the seeing of the worst input images, resulting in an even harsher trade-off between depth and seeing than for :ref:`Direct Coadds <direct_coadd>` and more correlated noise.  This is the only coadd for which nonlinear image combinations (such as a median or sigma-clipped mean) may be considered.

.. _stackfit:

StackFit
  An approach to source measurement (especially weak lensing shear estimation) that fits models to :ref:`Direct Coadds <direct_coadd>` after convolving with a PSF model generated using the :ref:`CoaddPsf <coaddpsf>` approach, developed by [Jee2011]_.  This avoids B-mode (and other) systematics that arise from poor modeling of PSF discontinuities, but is still lossy.  Contrast with :ref:`MultiFit <multifit>`.

.. _sufficient_statistic:

Sufficient Statistic
  Given a dataset and a likelihood that can be computed from it, a sufficient statistic for that dataset is any set of derived quantities from which the exact likelihood can also be computed.  In the context of this document, an optimal coadd is defined as any coadd that is a sufficient statistic for its input images for any likelihood that assumes a static (temporily nonvariable) sky.

.. _template:

Template
  A coadd used as the comparison image in difference imaging.  As the template must be convolved with a kernel that matches its PSF to that of the science image, :ref:`constant-PSF coadds <constant_psf_coadd>` are usually preferred, as they allow the matching kernel to be continuous.

.. _zackay_ofek_coadd:

Zackay/Ofek Coadd
  See :ref:`Kaiser Coadd <kaiser_coadd>`; from [Zackay2015]_, which indepenently derived Kaiser's result.


References
==========

.. [Szalay1999] `Szalay, Connolly, & Szokoly, 1999 <http://adsabs.harvard.edu/abs/1999AJ....117...68S>`_. *Simultaneous Multicolor Detection of Faint Galaxies in the Hubble Deep Field.* AJ, 117, 68.

.. [Jee2011] `Jee & Tyson, 2011 <http://adsabs.harvard.edu/abs/2011PASP..123..596J>`_. *Toward Precision LSST Weak-Lensing Measurement.* PASP, 123, 596.

.. [Kaiser2001] Kaiser, 2001.  *Addition of Images with Varying Seeing.* PSDC-002-011-xx.

.. [Zackay2015] `Zackay & Ofek, 2015 <http://adsabs.harvard.edu/abs/2015arXiv151206879Z>`_.  *How to coadd images? II. A coaddition image that is optimal for any purpose in the background dominated noise limit.* `arXiv:1512.06879 <http://arxiv.org/abs/1512.06879>`_
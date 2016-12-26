:tocdepth: 1


Notation and Conventions
========================

Given multiple overlapping images :math:`\{z_1({\bf r}), z_2({\bf r}), ...\}`, we designate by :math:`R({\bf r})` the set of image indices that overlap a point on the sky :math:`{\bf r}`.  We assume that all images have already been resampled to a common coordinate system and common photometric units.  Each image is related to the true sky :math:`f({\bf r})` via its point spread function (PSF) :math:`\phi_i({\bf r}, {\bf t})` by

.. math::
  z_i({\bf r}) = \int \! d^2{\bf t} \phi_i({\bf r}, {\bf t}) \, f({\bf t})
    + \mathrm{noise}
  :label: eq:psf_definition_continuous

It will be convenient to rewrite this integral as a sum, and hence define the PSF and sky as discrete quantities rather than continuous functions (note that :math:`{\bf r}` already only takes discrete values).  The true :math:`f` has power at arbitrarily high frequencies, so we cannot simply sample it on a grid.  But the PSF does not; we can always choose some grid with points :math:`{\bf s}` upon which we can sample the PSF and exactly reconstruct it with Sinc interpolation:

.. math::
  \phi_i({\bf r}, {\bf t}) = \sum_{\bf s} \phi_i({\bf r}, {\bf s})
    \, \mathrm{sinc}({\bf t} - {\bf s})

We can then insert this into :eq:`eq:psf_definition_continuous`, and reorder the sum and integral:

.. math::
  z_i({\bf r}) = \sum_{\bf t} \phi_i({\bf r}, {\bf t}) \int \! d^2{\bf s} \,
    \mathrm{sinc}({\bf t} - {\bf s}) \, f({\bf s})
    + \mathrm{noise}

This lets us identify the *sinc-convolved sky* :math:`h` as

.. math::
  h({\bf s}) \equiv \int\! d^2 {\bf r} \; \mathrm{sinc}({\bf r}-{\bf s}) \, f({\bf r})

Note that :math:`h` contains all in information about the true sky, and hence we can use it instead of :math:`f` to form a discrete analog of :eq:`eq:psf_definition_continuous`:

.. math::
  z_i({\bf r}) = \sum_s \phi_i({\bf r}, {\bf s}) \, h({\bf s})
    + \mathrm{noise}
  :label: eq:psf_definition

We will use this form for the remainder of the paper.

We assume all images have Gaussian noise described by a covariance matrix :math:`C_i({\bf r}, {\bf s})`.  Because astronomical images typically have approximately uncorrelated noise only in their original coordinate system, we should not assume that the covariance matrices in the common coordinate system are diagonal.

This notation uses function arguments for spatial indices and subscripts for quantities corresponding to different exposures, but (as described above) this does not imply that the spatial indices are continuous.  The spatial variables (typically :math:`{\bf r}` and :math:`{\bf s}`) should be assumed to take only discrete values, and indeed we will at times use matrix notation for sums over pixels (in which images are vectors, with the spatial index flattened):

.. math::
  {\bf z}_i = \left[
    \begin{array}{c}
      z_i(\{0,0\}) \\
      z_i(\{0,1\}) \\
      \vdots \\
      z_i(\{0,N\}) \\
      z_i(\{1,0\}) \\
      z_i(\{1,1\}) \\
      \vdots \\
      z_i(\{1,N\}) \\
      \vdots \\
      z_i(\{1,N\}) \\
      z_i(\{N,0\}) \\
      z_i(\{N,1\}) \\
      \vdots \\
      z_i(\{N,N\})
    \end{array}
  \right]

In matrix notation, :eq:`eq:psf_definition` is simply

.. math::
  {\bf z}_i = \boldsymbol{\phi}_i {\bf h} + \mathrm{noise}

Note that in matrix notation we continue to use subscripts to refer to exposure indices, not spatial indices; at no point will we use matrix notation to represent a sum over exposure indices.

We assume the true sky is the same in all input images (i.e. there is no variability).  This is clearly not true in practice, but coaddition algorithms are explicitly focused on capturing only the static sky.

Lossy Coaddition Algorithms
===========================

Direct Coadds
-------------

The simplest way to build a coadd is to simply form a linear combination of the images that overlap a point on the sky.  The coadd image :math:`z_\mathrm{dir}({\bf r})` is then

.. math::
  z_\mathrm{dir}({\bf r}) = \sum_{i \in R({\bf r})} w_i({\bf r}) \, z_i({\bf r})
  :label: eq:dir_coadd_def

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

In PSF-matched coaddition, input images are convolved by a kernel that matches their PSF to a predefined constant PSF before they are combined.  If :math:`\phi_\mathrm{pm}({\bf r})` is the predefined PSF for the coadd, then the matching kernel :math:`K_i({\bf r}, {\bf s})` is defined such that

.. math::
  \sum_{\bf u} \! K_i({\bf r}, {\bf u}) \, \phi_i({\bf u}, {\bf s})
    = \phi_\mathrm{pm}({\bf r}-{\bf s})

Typically :math:`K` is parametrized as a smoothly varying linear combination of basis functions.  The details of fitting it given a target coadd PSF and input image PSF models is beyond the scope of this document; see e.g. [Alard1998]_ for more information.

Because deconvolution is (at best) noisy, convolution with :math:`K_i` will generally increase the size of the PSF.  This highlights the big disadvantage of PSF-matched coadds: the images with the best seeing must be degraded to match a target PSF whose sizes is determined by the worst of the images to be included in the coadd.  Thus PSF-matched coadds must either include only the best-seeing images (sacrificing depth) or suffer from a worst-case coadd PSF.

After PSF-matching, the coadd is constructed in the same way as a direct coadd:

.. math::
  z_\mathrm{pm}({\bf r}) = \sum_{i \in R({\bf r})} w_i({\bf r}) \,
      \sum_{\bf u} K_i({\bf r}, {\bf u}) \, z_i({\bf u})
  :label: eq:pm_coadd_def

The PSF on the coadd is of course just :math:`\phi_\mathrm{pm}({\bf r})`, and the pixel covariance on the coadd is

.. math::
  C_\mathrm{pm}({\bf r}, {\bf s}) =
    \sum_{i \in \left[R({\bf r}) \cap R({\bf s})\right]} \!\!\!
        w_i({\bf r}) \, w_i({\bf s}) \,
        \sum_{\bf u} \sum_{\bf v} K_i({\bf r}, {\bf u}) \,
        K_i({\bf s}, {\bf v}) \,
        C_i({\bf u}, {\bf v})

Typically, the covariance terms in the uncertainty are simply ignored and only the variance is propagated, though this can result in a signficant misestimation of the uncertainty in measurements made on the coadd.


Outlier Rejection and Nonlinear Statistics
------------------------------------------

A common -- but often misguided -- practice in coaddition is to use a nonlinear statistic to combine pixels, substituting the weighted mean in :eq:`eq:dir_coadd_def` and :eq:`eq:pm_coadd_def` for a median or sigma-clipped mean.  The goal is to reject artifacts without explicitly detecting them on each image; the problem is that this assumes that the pixel values that go into a particular coadd pixel are drawn from distributions with the same mean.

This is not true when input images have different PSFs, as in direct coaddition.  Building a direct coadd with median or any amount of sigma-clipping will typically result in the cores of brighter stars being clipped in the the best seeing images, resulting in a flux-dependent (i.e. ill-defined) PSF.  Even extremely soft (e.g. 10-sigma) clipping is unsafe; the usual Gaussian logic concerning the number of expected outliers is simply invalid when the inputs are not drawn from the same distribution.

The presence of correlated noise means that even PSF-matched coadds cannot be built naively with nonlinear statistics.  In PSF-matched coadds, all pixels at the same point are drawn from distributions that have the same mean, but they are not identical distributions.  As a result, nonlinear statistics do not produce an ill-defined PSF when the inputs are PSF-matched, but their outlier rejection properties do not operate as one would naively expect, making it hard to predict how well any statistic will actually perform at eliminating artifacts (or not eliminating valid data).  Nonlinear statistics also make it impossible to correctly propagate uncertainty to coadds, as long as they are used to compute each coadd pixel independently.


Optimal Coaddition Algorithms
=============================

Likelihood Coadds
-----------------

An optimal coadd is one that is a :ref:`sufficient statistic <sufficient_statistic>` for the true sky: we can use it to compute the likelihood of a model of the true (static) sky, yielding the exact same computation as if we had computed the joint likelihood of that model over all the input images.  This joint likelihood is thus a natural starting point for deriving an optimal coadd.

The log likelihood of a single input image :math:`{\bf z}_i` is (in matrix notation)

.. math::
  L_i = -\frac{1}{2}
    \left[
      {\bf z}_i - \boldsymbol{\phi}_i{\bf h}
    \right]^T
    {\bf C}_i^{-1}
    \left[
      {\bf z}_i - \boldsymbol{\phi}_i{\bf h}
    \right]

The joint likelihood for all images is just the product of the per-image likelihoods since the images are independent. The joint log likelihood is thus the sum of the input log likelihoods:

.. math::
  L = \sum_i L_i
    = -\frac{1}{2} \sum_i \left[
          {\bf z}_i - \boldsymbol{\phi}_i{\bf h}
        \right]^T
        {\bf C}_i^{-1}
        \left[
          {\bf z}_i - \boldsymbol{\phi}_i{\bf h}
        \right]

By expanding this product, we can identify terms that include different powers of :math:`{\bf h}`:

.. math::
  L =& -\frac{1}{2} \sum_i {\bf z}_i^T {\bf C}_i^{-1} {\bf z}_i
    +  {\bf h}^T \left[
        \sum_i \boldsymbol{\phi}_i^T {\bf C}_i^{-1} {\bf z}_i
      \right]
    - \frac{1}{2} {\bf h}^T \left[
        \sum_i \boldsymbol{\phi}_i^T {\bf C}_i^{-1} \boldsymbol{\phi}_i
      \right] {\bf h} \\
    =& -\frac{k}{2}
    +  {\bf h}^T \boldsymbol{\Psi}
    - \frac{1}{2} {\bf h}^T \boldsymbol{\Phi} {\bf h}

with

.. math::
  k =& \sum_i {\bf z}_i^T {\bf C}_i^{-1} {\bf z}_i \\
  \boldsymbol{\Psi} =&
    \sum_i \boldsymbol{\phi}_i^T {\bf C}_i^{-1} {\bf z}_i \\
  \boldsymbol{\Phi} =&
    \sum_i \boldsymbol{\phi}_i^T {\bf C}_i^{-1} \boldsymbol{\phi}_i

These three terms represent a coadd of sorts.  :math:`\boldsymbol{\Psi}` is an image-like quantity, and :math:`\boldsymbol{\Phi}` behaves much like an (inverse) pixel covariance matrix.  Together with the scalar :math:`k` these are a sufficient statistic for :math:`{\bf h}`, and hence we can think of them as a form of optimal coadd, albeit one we cannot use in the usual way.  In particular, the covariance-like term :math:`\boldsymbol{\Phi}` does much more than just carry uncertainty information, as it captures what we typically think of as the PSF as well.  We will refer to the combination of :math:`\boldsymbol{\Psi}`, :math:`\boldsymbol{\Phi}`, and :math:`k` as a "likelihood coadd".

The fact that we cannot interpret a likelihood coadd in the same way as other astronomical images is inconvenient, but the real problem lies in its computational cost: :math:`\boldsymbol{\Phi}` is extremely large; while it is sparse, even just its nonzero elements would consume approximately 200GB in single precision for a single-patch 4k :math:`\times` 4k coadd.  While the same is broadly true of any detailed attempt to capture coadd uncertainty, :math:`\boldsymbol{\Phi}` has even more nonzero elements than :math:`{\bf C}_\mathrm{dir}` or :math:`{\bf C}_\mathrm{pm}`, and it plays a much more important role.  Approximating :math:`{\bf C}_\mathrm{dir}` and :math:`{\bf C}_\mathrm{pm}` generally implies incompletely or incorrectly propagating uncertainties, generally by a small amount, while approximating :math:`\boldsymbol{\Phi}` also implies incorrectly modeling the PSF.


Decorrelated Coadds
-------------------

The solution to the first problem of likelihood coadds -- that the images cannot be interpreted in the traditional way -- is to factor :math:`\boldsymbol{\Phi}`.  This is no small task, given the size of :math:`\boldsymbol{\Phi}`, but if it can be done, it also hints at a solution to the more serious computational problems with likelihood coadds.

Specifically, we assume a factorization of the form

.. math::
  \boldsymbol{\Phi} = \boldsymbol{\phi}_\mathrm{dec}^T
    {\bf C}_\mathrm{dec}^{-1}
    \boldsymbol{\phi}_\mathrm{dec}
  :label: eq:decorrelated_factorization

where :math:`\boldsymbol{\phi}_\mathrm{dec}` is a compact kernel and :math:`{\bf C}_\mathrm{dec}` is a nearly diagonal matrix.  Given that we have identified :math:`\boldsymbol{\Phi}` as representing the (inverse) covariance matrix of a likelihood coadd, this factorization essentially represents an attempt to *decorrelate* the noise on the likelihood coadd.  This is not quite sufficient, however; we also need to simultaneously solve for :math:`{\bf z}_\mathrm{dec}`` in

.. math::
  \boldsymbol{\Psi} = \boldsymbol{\phi}_\mathrm{dec}^T
    {\bf C}_\mathrm{dec}^{-1} {\bf z}_\mathrm{dec}
  :label: eq:decorrelated_coadd

As the notation implies, this allows us to rewrite the joint log likelihood as

.. math:: L = -\frac{1}{2}
  \left[
    {\bf z}_\mathrm{dec} - \boldsymbol{\phi}_\mathrm{dec} {\bf h}
  \right]^T
  {\bf C}_\mathrm{dec}^{-1}
  \left[
    {\bf z}_\mathrm{dec} - \boldsymbol{\phi}_\mathrm{dec} {\bf h}
  \right]

This identifies :math:`{\bf z}_\mathrm{dec}` as the decorrelated coadd image, :math:`{\bf C}_\mathrm{dec}` as its covariance matrix, and :math:`\boldsymbol{\phi}_\mathrm{dec}` as its PSF.  These can be used in exactly the same way as the corresponding single-exposure quantities.  As such, this is essentially the ideal coadd: it is formally optimal, can be used in the same way as any standard image, and makes no restrictive assumptions about the input images.

The problem is of course the computational expense.  Despite the fact that we have derived the decorrelated coadd from the likelihood coadd, we do not necessarily need to produce a full likelihood coadd first; it may be possible to devise an algorithm that factors :math:`\boldsymbol{\Phi}` in small regions as it is constructed.  And the decorrelated coadd quantities :math:`{\bf C}_\mathrm{dec}` and :math:`\boldsymbol{\phi}_\mathrm{dec}` may be much more amenable to compression than :math:`\boldsymbol{\Phi}`.

Because we have merely specified that :math:`{\bf C}_\mathrm{dec}` be "nearly" diagonal, this decomposition is not unique, and we have considerable flexibility to move power between :math:`\boldsymbol{\phi}_\mathrm{dec}` to make computation, storage, and use more efficient (without any change in the formal optimality).  Generally speaking, though, we want these quantities to mimic their standard image counterparts:

- We want :math:`{\bf C}_\mathrm{dec}` to be close to diagonal, and to capture small-scale changes in the variance due to bright objects.
- We want :math:`\boldsymbol{\phi}_\mathrm{dec}` to be compact and to vary smoothly (and slowly) over the image, to allow us to approximate the PSF as spatially constant over the scale of an object.

These are generally competing goals, as can be seen from the limiting cases (which are not necessarily solutions, especially when :eq:`eq:decorrelated_coadd` is considered)

.. math::
  {\bf C}_\mathrm{dec} =& \, \boldsymbol{\Phi}^{-1} \\
  \boldsymbol{\phi}_\mathrm{dec} =& \, {\bf I}

and

.. math::
  {\bf C}_\mathrm{dec} =& \, {\bf I} \\
  \boldsymbol{\phi}_\mathrm{dec} =& \, \boldsymbol{\Phi}^{1/2}

The former has a constant delta function PSF (recall that the pixel response is still embedded in the model) and highly correlated noise; the latter has white, uncorrelated noise and a non-compact PSF that can vary significantly from pixel to pixel.  Nevertheless, intuition suggests that it should be possible to achieve a solution in which the effective PSF is compact and fully continuous or piecewise continuous over large areas while the uncertainty is nondiagonal only in the neighborhood of boundary regions where the number of input images changes.

In addition to being familiar and hence convenient for downstream processing, optimizing these criteria should also make storage of :math:`{\bf C}_\mathrm{dec}` and :math:`\boldsymbol{\phi}_\mathrm{dec}` much more efficient.  Depending on how close to diagonal we can make it, :math:`{\bf C}_\mathrm{dec}` could require little more storage than the coadd image itself.  If we impose a continuous :math:`\boldsymbol{\phi}_\mathrm{dec}`, we can represent it as an interpolated function in essentially the same way we represent per-exposure PSF models.

Unfortunately, a general algorithm for computing the decorrelation factorization does not yet exist, making decorrelated coadds a mostly theoretical concept at present.  Some initial steps towards developing such an algorithm will be discussed in later sections.


Kaiser Coadds
-------------

If the input images to a likelihood coadd meet certain restrictive conditions, an algorithm developed by [Kaiser2001]_ (and rediscovered by [Zackay2015]_) can be used to build decorrelated coadd.  These conditions include:

- The noise in the input images must be white and uncorrelated.
- The PSFs of the input images must (individually) be spatially constant.
- The input images have no missing pixels, and the coadd area does not include any boundaries where the number of input images changes.

Under these conditions, :math:`\boldsymbol{\Phi}` has no spatial variation, giving it a particularly simple form in Fourier space:

.. math::
  \tilde{\Phi}({\bf u}, {\bf v})
    = \delta({\bf u}, {\bf v}) \sum_i \frac{
        \tilde{\phi}_i({\bf u}) \, \tilde{\phi}_i({\bf v})
      }{
        C_i
      }
    = \delta({\bf u}, {\bf v}) \sum_i \frac{
        \left|\tilde{\phi}_i({\bf u})\right|^2
      }{
        C_i
      }

(recall that :math:`C_i` is now just a scalar, as the variance is constant and there is no covariance).  Recognizing that the matrix products in :eq:`eq:decorrelated_factorization` are just convolutions when the products are spatially constant, the Fourier-space equivalent for Kaiser coadds is

.. math::
  \tilde{\Phi}({\bf u}, {\bf v}) =
    \tilde{\phi}_\mathrm{ksr}({\bf u})
    \,
    \left[ C_\mathrm{ksr}^{-1} \, \delta({\bf u}, {\bf v}) \right]
    \,
    \tilde{\phi}_\mathrm{ksr}({\bf v})

The solution is trivial (and unique, assuming a normalized PSF):

.. math::
  \tilde{\phi}_\mathrm{ksr}({\bf u}) =&
    \sqrt{
      \frac{
        \sum_i \left|\tilde{\phi}_i({\bf u})\right|^2 \, C_i^{-1}
      }{
        \sum_i C_i^{-1}
      }
    }\\
  C_\mathrm{ksr} =& \frac{1}{\sum_i C_i^{-1}}

The equivalent for :math:`\boldsymbol{\Psi}` and :eq:`eq:decorrelated_coadd` is

.. math::
  \tilde{\Psi}({\bf u}) = \sum_i
    \frac{
      \tilde{\phi}^*_i({\bf u}) \, \tilde{z}({\bf u})
    } {
      C_i
    }
  =
  \tilde{\phi}_\mathrm{ksr}({\bf u}) \,
    C_\mathrm{ksr}^{-1} \, \tilde{z}_\mathrm{ksr}({\bf u})

with solution

.. math::
  \tilde{z}_\mathrm{ksr}({\bf u}) =
    \frac{
      \sum_i \tilde{\phi}^*_i({\bf u}) \, \tilde{z}_i({\bf u}) \, C_i^{-1}
    }{
      \sqrt{
        \left[
          \sum_i \left| \tilde{\phi}_i({\bf u}) \right|^2 \, C_i^{-1}
        \right]
        \left[
          \sum_i C_i^{-1}
        \right]
      }
    }

This differs from [Zackay2015]_\'s Eqn. 7 because they have redefined the flux units of the coadd to achieve unit variance on the coadd.

The problem with the Kaiser algorithm is its assumptions, which are simply invalid for any realistic coadd.  While the noise in an input image may be white in the neighborhood of faint sources, most images contain brighter objects (and faint objects near brigher objects as well).  In addition, the noise is never uncorrelated once the image has been resampled to the coadd coordinate system.  The noise assumptions by themselves are not too restrictive, however; the Kaiser algorithm is not optimal when these conditions are not met, but we only care deeply about optimality in the neighborhood of faint sources.  And ignoring additional covariance due to warping is no different from our usual approach with direct coadds.

The assumptions that the PSFs and input image set are fixed are more problematic, but this still leaves room for the Kaiser algorithm to be used to build "per object" coadds, in which we build separate coadds each small region in the neighborhood of a single object, and reject any input image that do not fully cover that region.  This would likely necessitate coadding multiple regions multiple times (for overlapping objects), and it isn't as useful as a traditional coadd (especially considering that it can't be used for detection), but it may still have a role to play.

A more intriguing possibility is that the Kaiser approach could be used as one piece of a larger algorithm to build general decorrelated coadds.  One could imagine an iterative approach to solving :eq:`eq:decorrelated_factorization` and :eq:`eq:decorrelated_coadd` by minimizing a metric such as

.. math::
  q = \left|
        \boldsymbol{\Phi}
        - \boldsymbol{\phi}_\mathrm{dec}^T
          {\bf C}_\mathrm{dec}^{-1}
          \boldsymbol{\phi}_\mathrm{dec}
      \right|
    + \left|
        \boldsymbol{\Psi}
        - \boldsymbol{\phi}_\mathrm{dec}^T
          {\bf C}_\mathrm{dec}^{-1}
          {\bf z}_\mathrm{dec}
      \right|
    + \lambda \left|
        {\bf C}_\mathrm{dec}^{-1}
          - \mathrm{diag}({\bf C}_\mathrm{dec}^{-1})
      \right|

where :math:`\boldsymbol{\phi}_\mathrm{dec}` is parametrized as a smoothly-varying interpolation of a set of kernel basis functions, and :math:`\lambda` controls how strongly off-diagonal elements of :math:`{\bf C}_\mathrm{dec}^{-1}` are penalized.  This is a massive optimization problem if applied to a full coadd patch, but the structure of :math:`\boldsymbol{\Phi}` only indirectly couples pixels that are more than twice the PSF width apart; this suggests we could proceed by iteratively solving small regions independently -- if we have a good guess at an approximate solution.  The Kaiser algorithm provides exactly this: we can use the Kaiser method to estimate the PSF, and a diagonal covariance matrix at multiple points on the image, and then simply interpolate between them to generate our initial guess.  Just imposing the Kaiser PSF (or a small perturbation to it) as the final PSF may also be feasible.  This would only require us to solve for :math:`{\bf C}_\mathrm{dec}^{-1}` and :math:`{\bf z}_\mathrm{dec}`, dramatically reducing the scale of the problem.

Constant PSF Coadds
-------------------

A simple but potentially useful twist on the decorrelated coadd approach is to decorrelate only to a predefined constant PSF.  This would produce a coadd with many of the benefits of a PSF-matched coadd, but with no seeing restrictions on the input images and a much smaller final PSF.  Like a PSF-matched coadd, significant pixel correlations could remain in this scenario (it is unclear which approach would have more), but this coadd would enable the measurement of consistent colors and could also serve as a template for difference imaging.  Both of these are cases where having improved depth and a smaller PSF in the coadd could be critical.

Having a consistent PSF across bands is the only way to formally measure a consistent color, but using traditional PSF-matched coadds for this ensures these colors will have lower SNR than model-based measurements that operate on individual exposures (which are always at least somewhat biased).  If the constant-PSF coadd is instead generated using the decorrelated coadd approach, the SNR of consistent colors could be much more competitive.

The potential gains for difference imaging are even larger: the PSF size on the coadd puts a lower limit on the PSF size of an input exposure that can be differenced in it, which could require us to throw away or degrade our best images simply because we don't have a coadd good enough to difference with it. [#preconvolution]_  Difference imaging algorithms also become dramatically more complex when noise from the template cannot be neglected when compared with the noise in the exposure being differenced; this requires that the template have a large number of exposures.  This is challenging when traditional PSF-matched coaddition is used and the coadd PSF must be optimized along with the depth, and it may be even more challenging if mitigating chromatic PSF effects requires templates binned in airmass or some other approach that effectively adds new degrees of freedom to template generation.

.. [#preconvolution] The "preconvolution" approach to difference imaging decreases this lower limit (possibly to the point where it is unimportant), but is also an unproven technique.


Coadds for Source Detection
===========================

Detection Maps
--------------

The approach to source detection in LSST is derived from the likelihood of a single isolated point source of flux :math:`\alpha` centered on pixel :math:`\boldsymbol{\mu}`:

.. math::
  L =& -\frac{1}{2} \sum_i \sum_{{\bf r}, {\bf s}}
        \left[
          z_i({\bf r}) - \alpha\,\phi_i(\boldsymbol{\mu} - {\bf s})
        \right]
        \left[C_i^{-1}({\bf r}, {\bf s}) \right]
        \left[
          z_i({\bf s}) - \alpha\,\phi_i(\boldsymbol{\mu} - {\bf s})
        \right] \\
    =& -\frac{k}{2}
        + \alpha\Psi(\boldsymbol{\mu})
        - \frac{\alpha^2}{2}\Phi(\boldsymbol{\mu}, \boldsymbol{\mu})

At fixed :math:`\boldsymbol{\mu}`, we can solve for :math:`\alpha` by setting the first derivative of :math:`L` to zero:

.. math::
  \frac{\partial L}{\partial \alpha}
    = \Psi(\boldsymbol{\mu})
    - \alpha\Phi(\boldsymbol{\mu}, \boldsymbol{\mu})
    = 0

which yields

.. math::
  \hat{\alpha}(\boldsymbol{\mu})
    = \frac{
        \Psi(\boldsymbol{\mu})
      }{
        \Phi(\boldsymbol{\mu},\boldsymbol{\mu})
      }

Similarly, the variance in the flux can be computed from the inverse of the second derivative:

.. math::
  \sigma_{\alpha}^2(\boldsymbol{\mu})
    = \left( -\frac{\partial^2 L}{\partial \alpha^2} \right)^{-1}
    = \left[\Phi(\boldsymbol{\mu},\boldsymbol{\mu})\right]^{-1}

The point-source SNR at position :math:`\boldsymbol{\mu}` is then

.. math::
  \nu(\boldsymbol{\mu}) \equiv
  \frac{
    \hat{\alpha}(\boldsymbol{\mu})
  }{
    \sigma_{\alpha}(\boldsymbol{\mu})
  } =
    \frac{
        \Psi(\boldsymbol{\mu})
      }{
        \sqrt{\Phi(\boldsymbol{\mu},\boldsymbol{\mu})}
      }
  :label: eq:detection_map

To detect point sources, we simply threshold on :math:`\boldsymbol{\nu}`, which we call a *detection map*.   We can construct this from the components of a likelihood coadd with a crucial simplification: we only require the diagonal of :math:`\boldsymbol{\Phi}`, making what had been a computationally infeasible method quite practical.  This holds only because we have assumed an isolated point source, however; optimal detection of extended sources or blended sources would require at least some off-diagonal elements of :math:`\boldsymbol{\Phi}`.  In practice, we instead just look for multiple peaks in above-threshold regions in :math:`\boldsymbol{\nu}` as defined above, and bin the image to detect extended low-surface-brightness sources.


Optimal Multi-Band Detection
----------------------------

Just as optimal detection in monochromatic images requires that we know the signal of interest (a point source with a known PSF), optimal detection over multi-band observations requires that we know both the spectral energy distribution (SED) of the target objects and the bandpass.  More precisely, we need to know the integral of these quantities:

.. math::
  \beta_i = \int S(\lambda) \, T_i(\lambda) \, d\lambda

where :math:`T_i(\lambda)` is the normalized system response for observation :math:`i` and :math:`S(\lambda)` is the normalized SED of the target source.  The point source likelihood is then

.. math::
  L =& -\frac{1}{2} \sum_i \sum_{{\bf r}, {\bf s}}
        \left[
          z_i({\bf r}) - \alpha\,\beta_i\,\phi_i(\boldsymbol{\mu} - {\bf s})
        \right]
        \left[C_i^{-1}({\bf r}, {\bf s}) \right]
        \left[
          z_i({\bf s}) - \alpha\,\beta_i\,\phi_i(\boldsymbol{\mu} - {\bf s})
        \right] \\
    =& -\frac{k}{2}
        + \alpha\Psi_{\beta}(\boldsymbol{\mu})
        - \frac{\alpha^2}{2}\Phi_{\beta}(\boldsymbol{\mu}, \boldsymbol{\mu})

with

.. math::
  \boldsymbol{\Psi}_{\beta} =&
    \sum_i \beta_i \boldsymbol{\phi}_i^T {\bf C}_i^{-1} {\bf z}_i \\
  \boldsymbol{\Phi}_{\beta}=&
    \sum_i \beta_i^2 \boldsymbol{\phi}_i^T {\bf C}_i^{-1} \boldsymbol{\phi}_i

As the notation suggests, this is just a likelihood coadd with the inputs reweighted according to the target SED, and we can similarly form a detection map from it:

.. math::
  \nu_{\beta}(\boldsymbol{\mu}) =
      \frac{
          \Psi_{\beta}(\boldsymbol{\mu})
        }{
          \sqrt{\Phi_{\beta}(\boldsymbol{\mu},\boldsymbol{\mu})}
        }
  :label: eq:multiband_detection

In practice, the differences in throughput for different observations with the same bandpass is small enough to be neglected for detection purposes, and we could thus build :math:`\Phi_{\beta}` and :math:`\Psi_{\beta}` from per-band coadds of the standard :math:`\Phi` and :math:`\Psi`.  This makes it feasible to detect objects with unknown SEDs by quickly constructing detection maps for a library of proposed SEDs, and then merging those detections.

Chi-Squared Coadds
------------------

An alternate approach to multi-band coaddition developed by [Szalay1999]_ is to instead build a coadd that tests the null hypothesis that a pixel is pure sky.  While [Szalay1999]_ does not specify fully how to handle the spatial dimensions, we can combine their method with the likelihood coadd approach above.  This yields a detection map that is exactly the same as :eq:`eq:detection_map`, but with :math:`\Psi` and :math:`\Phi` summed over images from multiple bandpasses.  The probability distribution of :math:`\nu^2` is then a :math:`\chi^2` distribution, allowing the hypothesis test to be carried out by filtering on a monotonic function of the :math:`\nu`.

This is equivalent to setting :math:`\beta_i=1` in :eq:`eq:multiband_detection`, which is not the same as assuming a flat SED; in the background-dominated limit, it is actually the same as assuming that objects have the same SED as the sky.  From this perspective, it is clear that :math:`\chi^2` coadds are not formally optimal for the detection of most sources, but they may be close enough that detection on them with a slightly lower threshold may be more computationally efficient than trying a large library of proposed SEDs.


Quantitative Comparison
=======================

.. image:: /_static/comparison.png
   :target: ../../_static/comparison.png
   :alt: Comparison of seeing and depth for different coadd algorithms

The figure above shows predicted effective FWHM (calculated from PSF effective area) and 5-sigma point source limiting magnitude for different coadd algorithms.  Each data point in the histograms represents a coadd of 200 exposures, with seeing drawn from a log-normal distribution centered at 0.7 arcseconds and depth drawn from a normal distribution centered around 24.7.  Direct and PSF-matched coadds are weighted to optimize point source SNR.  Input PSFs use a Kolmogorov profile.

Note that the direct algorithm actually produces a smaller PSF than the Kaiser algorithm, even when the worst exposures are included (our choice of weight function strongly downweights these images).  This does *not* mean that it contains any more small-scale spatial information than the Kaiser coadd, as it always has lower SNR.  Even so, the improvement from direct to Kaiser algorithm is modest: when all exposures are included in the direct coadd, the Kaiser algorithm is only 0.1 magnitudes deeper.  The improvement from PSF-matched to direct coaddition is substantial in both PSF size and depth, especially when all exposures are included.  Imposing a cut on seeing percentile is clearly important for PSF-matched coaddition, but may not be important for direct coaddition, at least with the above choice of weight function.


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

.. [Alard1998] `Alard & Lupton, 1998 <http://adsabs.harvard.edu/abs/1998ApJ...503..325A>`_. *A Method for Optimal Image Subtraction.* ApJ, 503, 325.

.. [Szalay1999] `Szalay, Connolly, & Szokoly, 1999 <http://adsabs.harvard.edu/abs/1999AJ....117...68S>`_. *Simultaneous Multicolor Detection of Faint Galaxies in the Hubble Deep Field.* AJ, 117, 68.

.. [Jee2011] `Jee & Tyson, 2011 <http://adsabs.harvard.edu/abs/2011PASP..123..596J>`_. *Toward Precision LSST Weak-Lensing Measurement.* PASP, 123, 596.

.. [Kaiser2001] Kaiser, 2001.  *Addition of Images with Varying Seeing.* PSDC-002-011-xx.

.. [Zackay2015] `Zackay & Ofek, 2015 <http://adsabs.harvard.edu/abs/2015arXiv151206879Z>`_.  *How to coadd images? II. A coaddition image that is optimal for any purpose in the background dominated noise limit.* `arXiv:1512.06879 <http://arxiv.org/abs/1512.06879>`_

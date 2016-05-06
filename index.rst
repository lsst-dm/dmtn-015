..
  Content of technical report.

  See http://docs.lsst.codes/en/latest/development/docs/rst_styleguide.html
  for a guide to reStructuredText writing.

  Do not put the title, authors or other metadata in this document;
  those are automatically added.

  Use the following syntax for sections:

  Sections
  ========

  and

  Subsections
  -----------

  and

  Subsubsections
  ^^^^^^^^^^^^^^

  To add images, add the image file (png, svg or jpeg preferred) to the
  _static/ directory. The reST syntax for adding the image is

  .. figure:: /_static/filename.ext
     :name: fig-label
     :target: http://target.link/url

     Caption text.

   Run: ``make html`` and ``open _build/html/index.html`` to preview your work.
   See the README at https://github.com/lsst-sqre/lsst-report-bootstrap or
   this repo's README for more info.

   Feel free to delete this instructional comment.

:tocdepth: 2

Lossy Coaddition Algorithms
===========================

Exact Coaddition Algorithms
===========================

Likelihood Coadds
-----------------

Decorrelated Coadds
-------------------

Kaiser Coadds
-------------

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
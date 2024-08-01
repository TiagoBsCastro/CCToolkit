.. CCToolkit documentation master file, created by
   sphinx-quickstart on Thu Aug  1 18:23:49 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CCToolkit Documentation
=======================

.. contents:: Table of Contents
   :depth: 2
   :local:

CCToolkit is a Python package designed for Cluster Cosmology calculations, including handling cosmological parameters, power spectrum computations, halo mass functions, and halo bias. It integrates with the CAMB library to deliver precise cosmological data processing.

Features
========

- **Cosmological Calculations**: Easily compute various cosmological quantities, including background quantities and power spectra.
- **Halo Mass Function (HMF)**: Calculate the multiplicity function and HMF parameters for different halo finders based on the model presented in `Castro et al. 2022 <https://inspirehep.net/literature/2132031>`_.
- **Halo Bias**: Includes functions to compute the linear halo bias, with corrections based on the Peak Background Split (PBS) model presented in Castro et al. 2024.
- **Utility Functions**: Provides useful utilities for the manipulation of the power spectrum.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   examples
   api/modules

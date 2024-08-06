Examples
========

To use CCToolkit, start by importing the necessary modules and initializing the ``CosmologyCalculator`` with your desired cosmological parameters:

.. code-block:: python

   from cctoolkit import CosmologyCalculator

   # Define cosmological parameters
   params = {
       'H0': 70.0,                # Hubble parameter at z=0 in km/s/Mpc
       'Ob0': 0.05,               # Baryon density parameter
       'Om0': 0.3,                # Total matter density parameter (Ob0 + Om_c0)
       'sigma8': 0.8,             # rms density fluctuation amplitude at 8 h^-1 Mpc
       'ns': 0.96,                # Scalar spectral index
       'TCMB': 2.7255,            # CMB Temperature at z=0 in K
       'mnu': 0.06,               # Sum of neutrino masses (eV)
       'num_massive_neutrinos': 1 # Number of massive neutrinos species
   }

   # Initialize the calculator
   cosmo_calc = CosmologyCalculator(params)

Halo Mass Function
------------------

Calculate the halo mass function for a range of masses:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   masses = np.logspace(13, 15.5, num=100)
   hmf = cosmo_calc.dndlnM(masses, 0)
   plt.loglog(masses, hmf)
   plt.xlabel(r"$M_{\rm vir}\,[M_\odot h^{-1}]$")
   plt.ylabel(r"$\frac{{\rm d} n}{{\rm d} \log M}\,[{\rm Mpc}^{-3} h^{3}]$")
   plt.show()

Halo Bias
---------

Compute the linear halo bias or the at redshift z = 0 for a given mass. CCToolkit can compute both the PBS prescription as well as the corrected model following Castro et al. 2024.

.. code-block:: python

   pbs = cosmo_calc.pbs_bias(masses, 0)
   bias = cosmo_calc.bias(masses, 0)
   plt.loglog(masses, pbs, label=r'${\rm PBS}$')
   plt.loglog(masses, bias, label=r'${\rm Castro\, et\,al.\,2022}$')
   plt.xlabel(r"$M_{\rm vir}\,[M_\odot h^{-1}]$")
   plt.ylabel(r"$b(M)$")
   plt.legend()
   plt.show()


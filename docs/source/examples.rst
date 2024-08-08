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

Calculate the `Castro et al. 2023 <https://inspirehep.net/literature/2132031>`__ halo mass function for a range of masses:

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

Compute the linear halo bias or the at redshift z = 0 for a given mass. CCToolkit can compute both the PBS prescription as well as the corrected model following Castro et al. in prep.

.. code-block:: python

   pbs = cosmo_calc.pbs_bias(masses, 0)
   bias = cosmo_calc.bias(masses, 0)
   plt.loglog(masses, pbs, label=r'${\rm PBS}$')
   plt.loglog(masses, bias, label=r'${\rm Castro\, et\,al.\,2022}$')
   plt.xlabel(r"$M_{\rm vir}\,[M_\odot h^{-1}]$")
   plt.ylabel(r"$b(M)$")
   plt.legend()
   plt.show()

Impact of Baryons on Cluster and Group masses
---------------------------------------------

Compute the dark-matter only equivalent of a hydro-dynamic simulated object. CCToolkit can compute both the dark-matter only mass at the threshold of equivalence according to the quasi-adibatic model presented in `Castro et al. 2024 <https://inspirehep.net/literature/2718844>`__ but also to convert this to the virial definition assuming a concentration mass relation. For any realistic feedback mechanism, the last step is not very sensitive to the assumed mass-concentration relation. The assumed mass-concentration relation is given by `Duffy et al. 2008 <https://inspirehep.net/literature/783522>`__.

.. code-block:: python

    from cctoolkit import baryons
    from cctoolkit.cosmology import CosmologyCalculator
    z = 0.0
    # Magneticum cosmology
    params = {"Om0": 0.272, "Ob0": 0.272 * 0.168, "H0": 70.4, "ns": 0.963, "mnu":0, "num_massive_neutrinos": 0, "sigma8": 0.809}
    cosmo_calc = CosmologyCalculator()
    # Array of virial masses on hydro
    M = np.geomspace(1e13, 3e14)
    # Calculating the equivalent DMO mass at the threshold of equivalence according to the quasi-adiabatic model
    M_Delta_dmo, Delta = baryons.compute_dmo_mass(M, z, 0.168)
    Delta *= cctoolkit.utils.virial_Delta(cosmo_calc.Omega_m(z))
    # Converting to the virial mass
    mdmo = [baryons.compute_rec_mass(cosmo_calc, m, d, z) for m, d in zip(M_Delta_dmo, Delta)]

Using a tabulated matter power-spectrum
---------------------------------------

The ``CosmologyCalculator`` can also receive a tabulated power-spectrum. This is useful when analysing simulations which, due to backscaling, might have a realized power-spectrum that is not compatible with camb.

.. code-block:: python

    from cctoolkit.cosmology import CosmologyCalculator
    # Calling colossus to produce a tabulated Pk
    from colossus.cosmology import cosmology

    params = {'flat': True, 'H0': 67.321, 'Om0': 0.3158, 'Ob0': 0.0494, 'sigma8': 0.8102, 'ns': 0.9661}
    cosmo = cosmology.setCosmology("C0", params)
    params = {
       'H0': cosmo.H0,           # Hubble parameter at z=0 in km/s/Mpc
       'Ob0': cosmo.Ob0,         # Physical baryon density parameter
       'Om0': cosmo.Om0,         # Physical total matter density parameter (Ob0 + Om_c0)
       'sigma8': cosmo.sigma8,   # rms density fluctuation amplitude at 8 h^-1 Mpc
       'ns': cosmo.ns,           # Scalar spectral index
       'mnu': 0.0,               # Sum of neutrino masses (eV)
       'num_massive_neutrinos': 0,
    }
    M = np.geomspace(1e13, 1e16, 200)
    k = np.geomspace(1e-3, 1e1, 200)
    Pk = cosmo.matterPowerSpectrum(k, 0)
    cosmo_calc = CosmologyCalculator(params, power_spectrum=[k, Pk])
    dndlnM = cosmo_calc.dndlnM(M, 0)

Notice that simulations frequently ignores the radiation contribution. As we use camb as out backend for the ``CosmologyCalculator``, we can not produce a background without radiation. However, we can make its contribution insignificant setting the CMB temperature today to an unrealistic low value.

.. code-block:: python

    from cctoolkit.cosmology import CosmologyCalculator
  
    params = {'flat': True, 'H0': 67.321, 'Om0': 0.3158, 'Ob0': 0.0494, 'sigma8': 0.8102, 'ns': 0.9661, 'mnu': 0, 'num_massive_neutrinos': 0}
    cosmo = cosmology.setCosmology("C0", params)
    params['mnu'] = 0
    params['num_massive_neutrinos'] = 0
    params['TCMB'] = 0.5
    cosmo_calc = CosmologyCalculator(params)

When using a tabulated power-spectrum, CosmoCalculator will compute the growth factor solving Eq. (11) of `Linder and Jenkins 2003 <https://inspirehep.net/literature/618898>`__. If TCMB is lower than unity, EdS initial conditions are assumed at high-redshift. Otherwise, a radiation dominated solution is assumed.
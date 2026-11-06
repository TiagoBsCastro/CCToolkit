Examples
========

To use CCToolkit, start by importing the necessary modules and initializing the ``CosmologyCalculator`` with your desired cosmological parameters:

.. code-block:: python

   from cctoolkit import CosmologyCalculator

   # Define cosmological parameters
   params = {
       'H0': 70.0,                 # Hubble parameter at z=0 in km/s/Mpc
       'Ob0': 0.05,                # Baryon density parameter
       'Om0': 0.3,                 # Total matter density parameter (Ob0 + Om_c0)
       'sigma8': 0.8,              # rms density fluctuation amplitude at 8 h^-1 Mpc
       'ns': 0.96,                 # Scalar spectral index
       'TCMB': 2.7255,             # CMB Temperature at z=0 in K
       'mnu': 0.06,                # Sum of neutrino masses (eV)
       'num_massive_neutrinos': 1, # Number of massive neutrinos species
       'w0' : -1.1,                # Dark energy equation of state value today
       'wa' : 0.8                  # Dark energy equation of state rate wrt to a
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

Alternatively, one can specify the model `Castro et al. 2025 <https://inspirehep.net/literature/2910564>`__ for cosmologies with dynamic dark energy with:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   masses = np.logspace(13, 15.5, num=100)
   hmf = cosmo_calc.dndlnM(masses, 0, model='castro25')
   plt.loglog(masses, hmf)
   plt.xlabel(r"$M_{\rm vir}\,[M_\odot h^{-1}]$")
   plt.ylabel(r"$\frac{{\rm d} n}{{\rm d} \log M}\,[{\rm Mpc}^{-3} h^{3}]$")
   plt.show()

Halo Bias
---------

Compute the linear halo bias or the at redshift z = 0 for a given mass. CCToolkit can compute both the PBS prescription as well as the corrected model following `Castro et al. 2024b <https://inspirehep.net/literature/2824892>`__

.. code-block:: python

   pbs = cosmo_calc.pbs_bias(masses, 0)
   bias = cosmo_calc.bias(masses, 0)
   plt.loglog(masses, pbs, label=r'${\rm PBS}$')
   plt.loglog(masses, bias, label=r'${\rm Castro\, et\,al.\,2022}$')
   plt.xlabel(r"$M_{\rm vir}\,[M_\odot h^{-1}]$")
   plt.ylabel(r"$b(M)$")
   plt.legend()
   plt.show()

Nonlinear Halo Bias
-------------------

CCToolkit implements the DUCA III model (Castro 2026) of the nonlinear halo-matter bias :math:`b_{\rm hm}(k) = P_{\rm hm}(k)/P_{\rm mm}(k)` and of the halo auto power spectrum :math:`P_{\rm hh}(k)`, including stochasticity. It is calibrated on the DUCA simulations of `Castro et al. 2025 <https://inspirehep.net/literature/2910564>`__ for :math:`z < 2.5`, in :math:`w_0w_a` CDM cosmologies with massive neutrinos, up to :math:`k = 1.13`-:math:`2.25\,h\,{\rm Mpc}^{-1}` depending on the halo mass. On large scales :math:`b_{\rm hm} \to b_1`, the linear bias of `Castro et al. 2024b <https://inspirehep.net/literature/2824892>`__ recalibrated on the same simulations. All spectra refer to the cold matter (CDM plus baryons), the default ``var='cdm+b'`` of the calculator.

Up to :math:`k = 0.4\,h\,{\rm Mpc}^{-1}` the cross bias needs only the cosmology. The auto spectrum, and the cross bias on smaller scales, also need the nonlinear cold-matter power spectrum :math:`P_{\rm mm}(k)`, for instance from an emulator or, as below, from CAMB:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import camb
   from cctoolkit import CosmologyCalculator

   params = {'H0': 67.0, 'Ob0': 0.049, 'Om0': 0.315, 'As': 2.1e-9, 'ns': 0.965,
             'mnu': 0.06, 'num_massive_neutrinos': 1, 'w0': -0.9, 'wa': -0.3}
   cosmo_calc = CosmologyCalculator(params, zmax=2.5)
   z = 0.5
   k = np.geomspace(1e-3, 1.5, 300)

   # Cross bias of halos of given masses on large scales: shape (n_mass, n_k)
   masses = np.array([2e12, 1e13, 1e14, 1e15])
   b_hm, _ = cosmo_calc.nonlinear_bias(k[k <= 0.4], masses, z)

   # Nonlinear cold-matter (no neutrinos) power spectrum, here from CAMB's HMcode
   h = params['H0'] / 100
   pars = camb.set_params(H0=params['H0'], ombh2=params['Ob0'] * h**2,
                          omch2=(params['Om0'] - params['Ob0']) * h**2, mnu=params['mnu'],
                          num_massive_neutrinos=1, w=params['w0'], wa=params['wa'],
                          dark_energy_model='ppf', As=params['As'], ns=params['ns'],
                          redshifts=[z], kmax=10, NonLinear=camb.model.NonLinear_both)
   P_nl = camb.get_results(pars).get_matter_power_interpolator(
       nonlinear=True, var1='delta_nonu', var2='delta_nonu', hubble_units=True, k_hunit=True)
   P_mm = P_nl.P(z, k)

   # Halos in narrow mass bins (as in the calibration): b_hm and P_hh, shape (n_bins, n_k).
   # The number density of each bin follows from the Castro et al. 2025 mass function.
   edges = np.exp(np.arange(np.log(2e12), np.log(2e15), 0.2))
   b_hm, P_hh = cosmo_calc.nonlinear_bias(k, edges, z, matter_power=P_mm, bins=True)

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
   for i in range(0, len(edges) - 1, 10):
       ax1.semilogx(k, b_hm[i])
       ax2.loglog(k, P_hh[i])
   ax1.set_xlabel(r"$k\,[h\,{\rm Mpc}^{-1}]$")
   ax1.set_ylabel(r"$b_{\rm hm}(k)$")
   ax2.set_xlabel(r"$k\,[h\,{\rm Mpc}^{-1}]$")
   ax2.set_ylabel(r"$P_{\rm hh}(k)\,[h^{-3}\,{\rm Mpc}^3]$")
   plt.show()

``nonlinear_bias`` returns arrays of shape ``(n_mass, n_k)``; ``P_hh`` is ``None`` when the matter power spectrum or the number density is missing. With ``bins=True`` the masses are bin edges: the bias is averaged over each bin with the Castro et al. 2025 mass function, which also gives the number density. Pass ``number_density`` to use the abundance of your own halo sample, and ``return_variables=True`` to get the halo and cosmology descriptors used by the model (``halo_descriptors``). ``cctoolkit.NonlinearHaloBias`` evaluates the model directly on these descriptors, without CAMB:

.. code-block:: python

   from cctoolkit import NonlinearHaloBias

   model = NonlinearHaloBias()
   b_hm, P_hh, d = cosmo_calc.nonlinear_bias(k, edges, z, matter_power=P_mm, bins=True,
                                             return_variables=True)
   b1 = model.linear_bias(d['peak_height'], d['sigma_slope'], d['omega_cb_z'], d['S8'],
                          d['pbs_bias'])

Domain of validity:

- :math:`z < 2.5`, halos above 150 particles of the DUCA runs (for the Fiducial cosmology, :math:`M \gtrsim 1.5\times10^{12}\,M_\odot h^{-1}` from the high-resolution runs), and the cosmological parameter ranges of the DUCA I simulations.
- The model describes a halo sample by one set of properties, and it is calibrated on bins of width :math:`\Delta\ln M = 0.2`. A wider bin is therefore never evaluated at a single mass: ``nonlinear_bias`` splits it into sub-bins of at most that width and averages their :math:`b_{\rm hm}` with number-density weights, which is exact for :math:`P_{\rm hm}`. The auto spectrum of a wide bin would also need the cross spectra between its sub-bins, so :math:`P_{\rm hh}` is only returned for bins of at most :math:`\Delta\ln M = 0.2`.
- In these bins the halo shot noise dominates :math:`P_{\rm hh}` at :math:`k \gtrsim 0.3\,h\,{\rm Mpc}^{-1}`. For much denser samples (for instance a ``number_density`` far above the mass-function value), :math:`P_{\rm hh}` on small scales is an extrapolation.
- Between 0.4 and :math:`0.6\,h\,{\rm Mpc}^{-1}` the large-scale model is joined to its small-scale extension. Validate the join before transforming the predictions to configuration space.

Impact of Baryons on Cluster and Group masses
---------------------------------------------

Compute the dark-matter only equivalent of a hydro-dynamic simulated object. CCToolkit can compute both the dark-matter only mass at the threshold of equivalence according to the quasi-adibatic model presented in `Castro et al. 2024 <https://inspirehep.net/literature/2718844>`__ but also to convert this to the virial definition assuming a concentration mass relation. For any realistic feedback mechanism, the last step is not very sensitive to the assumed mass-concentration relation. The assumed mass-concentration relation is given by `Duffy et al. 2008 <https://inspirehep.net/literature/783522>`__.

.. code-block:: python

    import cctoolkit
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

Notice that simulations frequently ignores the radiation contribution. As we use camb as our backend for the ``CosmologyCalculator``, we can not produce a background without radiation. However, we can make its contribution insignificant setting the CMB temperature today to an unrealistic low value.

.. code-block:: python

    from cctoolkit.cosmology import CosmologyCalculator
  
    params = {'flat': True, 'H0': 67.321, 'Om0': 0.3158, 'Ob0': 0.0494, 'sigma8': 0.8102, 'ns': 0.9661, 'mnu': 0, 'num_massive_neutrinos': 0, 'TCMB': 0.5}
    cosmo_calc = CosmologyCalculator(params)

When using a tabulated power-spectrum, CosmoCalculator will compute the growth factor solving Eq. (11) of `Linder and Jenkins 2003 <https://inspirehep.net/literature/618898>`__. If TCMB is lower than unity, EdS initial conditions are assumed at high-redshift. Otherwise, a radiation dominated solution is assumed.

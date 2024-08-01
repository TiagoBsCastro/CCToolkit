"""
cosmology.py

This module defines the CosmologyCalculator class, which manages cosmological
parameters and calculations related to the power spectrum and other cosmological
quantities using CAMB.
"""
import numpy as np
try:
    from camb import CAMBparams, get_results, model
except ImportError:
    CAMBparams = None
    get_results = None
    model = None
from .utils import compute_sigma8_norm, ssq_given_Dk, d_s_given_Dk
from .hmf import multiplicity_function, best_fit_values_AHF, best_fit_values_ROCKSTAR, best_fit_values_SUBFIND, best_fit_values_VELOCIraptor
from .bias import corrected_bias

best_fits = {'AHF': best_fit_values_AHF,
             'ROCKSTAR': best_fit_values_ROCKSTAR,
             'SUBFIND': best_fit_values_SUBFIND,
             'VELOCIraptor': best_fit_values_VELOCIraptor}

class CosmologyCalculator:
    """
    A class to handle cosmological calculations using CAMB.

    This class facilitates the calculation of various cosmological quantities, including
    power spectra, density parameters, and functions related to the Halo Mass Function (HMF).
    It includes methods for calculating the Peak Background Split (PBS) bias and corrected bias,
    utilizing the CAMB library.

    Attributes:
    -----------
    params : dict
        Cosmological parameters for the model, including 'H0', 'Ob0', 'Om0', 'sigma8', 'ns', 'As', 'tau', 'TCMB', and 'mnu'.
    cosmo : object
        An instance of CAMB results, containing cosmological data.
    k : np.ndarray
        Array of wavenumbers used in the power spectrum calculations.
    Pk : dict
        Dictionary where keys are redshifts and values are arrays of the power spectrum values.
    Dk : dict
        Dictionary where keys are redshifts and values are arrays of dimensionless power spectrum values.
    z_vals : np.ndarray
        Array of redshift values used for power spectrum calculations.

    Methods:
    --------
    delta_c(Om, z=None):
        Calculates the critical density contrast (delta_c) for a given matter density.
    get_power_spectrum(z):
        Returns the wavenumbers and dimensionless power spectrum at a given redshift.
    set_cosmology(new_params):
        Updates the cosmological parameters and recalculates the power spectrum normalization.
    Omega_m(z):
        Returns the matter density parameter Omega_m as a function of redshift.
    Omega_nu(z):
        Returns the neutrino density parameter Omega_nu as a function of redshift.
    Omega_k(z):
        Returns the curvature density parameter Omega_k as a function of redshift.
    Omega_DE(z):
        Returns the dark energy density parameter Omega_DE as a function of redshift.
    critical_density(z):
        Returns the critical density at a given redshift.
    H(z):
        Returns the Hubble parameter H at a given redshift.
    peak_height(M, z):
        Calculates the peak-height delta_c / sigma(R(M), z) for a given mass and redshift.
    lagrangian_radius(M):
        Calculates the radius of the Lagrangian patch corresponding to a given mass M.
    sigma(R, z):
        Calculates the RMS fluctuation sigma(R, z) for a given radius and redshift.
    dlnsigma_dlnR(R):
        Calculates the derivative of the logarithm of sigma with respect to the logarithm of R.
    vfv(M, z, return_variables=False, halo_finder=best_fit_values_ROCKSTAR):
        Calculates the multiplicity function νf(ν) and related variables for a given mass and redshift.
    dndlnM(M, z, halo_finder='ROCKSTAR'):
        Calculates the halo mass function dn/dlnM, representing the number density of halos per logarithmic mass interval.
    pbs_bias(M, z, halo_finder='ROCKSTAR', return_variables=False):
        Calculates the PBS bias for halos of mass M at redshift z.
    bias(M, z, halo_finder='ROCKSTAR'):
        Calculates the corrected linear halo bias for given mass and redshift, using the PBS prediction and correction factors.

    Notes:
    ------
    - This class requires CAMB to be installed and available for the power spectrum calculations.
    - The methods assume that the input masses (M) are in units of solar masses over h (M_sun/h), and distances (R) are in Mpc/h.
    - The halo finder parameters must be provided correctly to ensure accurate mass function and bias calculations.
    """

    def __init__(self, params, power_spectrum=None, zmax=2.0, nz=100, var='cdm+b'):
        """
        Initialize the CosmologyCalculator with given cosmological parameters.

        Parameters:
        -----------
        params : dict
            Dictionary containing the cosmological parameters, including 'H0', 'Ob0', 'Om0', 'sigma8',
            'ns', 'As', 'tau', 'TCMB', and 'mnu'.
        power_spectrum : array, optional
            Array containing the power spectrum data (k and Pk). If None, the power spectrum
            will be computed from the parameters using CAMB.
        zmax : float, optional
            Maximum redshift for the calculation. Default is 2.0.
        nz : int, optional
            Number of redshift points. Default is 100.
        var : str, optional
            Which variable to use to compute the matter power spectrum. Default is 'cdm+b' (no neutrinos).
        """
        self.params = params
        self.zmax = zmax
        self.nz = nz
        self._z_vals_inv = np.linspace(0, zmax, nz)[::-1]

        if var == 'cdm+b':
            self._transfer_var = 8
        elif var == 'tot':
            self._transfer_var = 7
        else:
            raise RuntimeError(f"Not supported matter power-spectrum for var={var}!")

        if power_spectrum is not None:
            if (self.zmax != 0) or (self.nz != 1):
                raise RuntimeError("When using a custom power-spectrum zmax should be 0 and nz 1!\n \
                                   Control the redshift evolution changing sigma8 accordingly.")
            self._load_power_spectrum(power_spectrum)
            self._initialize_camb(params, background_only=True)
        else:
            self._initialize_camb(params)

    def _initialize_camb(self, params, background_only=False):
        """
        Initialize CAMB with the given cosmological parameters.

        Parameters:
        -----------
        params : dict
            Dictionary containing the cosmological parameters.
        background_only: bool
            Whether matter power-spectrum is expected to be computed or not.
        """
        if CAMBparams is None or get_results is None:
            raise ImportError("CAMB is not installed or not available.")
        
        camb_params = CAMBparams()
        camb_params.set_cosmology(
            H0=params['H0'],
            ombh2=params['Ob0'] * (params['H0'] / 100) ** 2,
            omch2=(params['Om0'] - params['Ob0']) * (params['H0'] / 100) ** 2,
            TCMB=params.get('TCMB', 2.7255),  # Default CMB temperature is 2.7255 K
            mnu=params.get('mnu', 0.06),  # Default sum of neutrino masses in eV
            num_massive_neutrinos=params.get('num_massive_neutrinos', 1)  # Number of massive neutrinos
        )
        camb_params.set_dark_energy(w=params.get('w0', -1.0), 
                                    cs2=1.0, 
                                    wa=params.get('wa', 0.0), dark_energy_model='ppf')
        if not background_only:
            camb_params.InitPower.set_params(As=params.get('As', 2e-9), ns=params['ns'])
            camb_params.set_matter_power(redshifts=self._z_vals_inv, kmax=10.0, k_per_logint=10)
        results = get_results(camb_params)
        self._cosmo = results

        if not background_only:
            self.k, self.z_vals, Pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=200, var1=self._transfer_var, var2=self._transfer_var)
            self._Dk = Pk * self.k ** 3 / (2*np.pi**2)
            self._Dk /= compute_sigma8_norm(self.k, self._Dk[0], self.params['sigma8'])
            self.Dk = {z: self._Dk[i] for i, z in enumerate(self.z_vals)}

    def _load_power_spectrum(self, power_spectrum):
        """
        Load the power spectrum data.

        Parameters:
        -----------
        power_spectrum_file : str
            Path to the file containing the power spectrum data (k and Pk).
        """
        self.k, Pk = power_spectrum
        self._Dk = Pk * self.k ** 3 / (2 * np.pi**2)
        self._Dk /= compute_sigma8_norm(self.k, self._Dk, self.params['sigma8'])
        self.Dk = {0: self._Dk}

    def delta_c(self, Om, z=None):
        """
        Calculate the critical density contrast (delta_c) using Bryan and Norman fit.

        Parameters:
        -----------
        Om : float or callable
            Matter density parameter. If callable, it should return Om as a function of z.
        z : float or array-like, optional
            Redshift at which to evaluate Om if Om is a function.

        Returns:
        --------
        float
            The critical density contrast delta_c.
        """
        if z is None:
            return 3.0 / 20.0 * (12.0 * np.pi) ** (2.0 / 3.0) * (1.0 + 0.012299 * np.log10(Om))
        else:
            return 3.0 / 20.0 * (12.0 * np.pi) ** (2.0 / 3.0) * (1.0 + 0.012299 * np.log10(Om(z)))
    
    def get_power_spectrum(self, z=0):
        """
        Get the wavenumber and dimensionless power spectrum arrays at a given redshift.

        Parameters:
        -----------
        z : float, optional
            Redshift at which to get the power spectrum. Default is 0.

        Returns:
        --------
        tuple
            (k, Dk) where k is the array of wavenumbers and Dk is the dimensionless power spectrum.
        """
        if z in self.Dk:
            return self.k, self.Dk[z]
        else:
            raise ValueError(f"Power spectrum for redshift z={z} not available.")

    def set_cosmology(self, new_params):
        """
        Update cosmological parameters and recalculate the power spectrum normalization.

        Parameters:
        -----------
        new_params : dict
            Dictionary of new cosmological parameters.
        """
        self.params = new_params
        self._initialize_camb(new_params)

    def Omega_m(self, z):
        """
        Get the matter density parameter Omega_m at a given redshift.

        Parameters:
        -----------
        z : float
            Redshift.

        Returns:
        --------
        float
            The matter density parameter Omega_m at redshift z.
        """
        return self._cosmo.get_Omega('cdm', z) + self._cosmo.get_Omega('baryon', z)
    
    def Omega_nu(self, z):
        """
        Get the neutrino density parameter Omega_nu at a given redshift.

        Parameters:
        -----------
        z : float
            Redshift.

        Returns:
        --------
        float
            The neutrino density parameter Omega_nu at redshift z.
        """
        return self._cosmo.get_Omega('nu', z)
    
    def Omega_k(self, z):
        """
        Get the curvature density parameter Omega_k at a given redshift.

        Parameters:
        -----------
        z : float
            Redshift.

        Returns:
        --------
        float
            The curvature density parameter Omega_k at redshift z.
        """
        return self._cosmo.get_Omega('k', z)
    
    def Omega_DE(self, z):
        """
        Get the dark energy density parameter Omega_DE at a given redshift.

        Parameters:
        -----------
        z : float
            Redshift.

        Returns:
        --------
        float
            The dark energy density parameter Omega_DE at redshift z.
        """
        return self._cosmo.get_Omega('de', z)

    def critical_density(self, z):
        """
        Get the critical density at a given redshift.

        Parameters:
        -----------
        z : float
            Redshift.

        Returns:
        --------
        float
            The critical density at redshift z in units of 10^10 M_sun / (h^2 Mpc^3).
        """
        H0 = self.params['H0']  # Hubble parameter today in km/s/Mpc
        H_z = self._cosmo.hubble_parameter(z)  # H(z) in km/s/Mpc
        G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
        Mpc_to_m = 3.085677581491367e22  # Conversion from Mpc to meters
        Msun_to_kg = 1.98847e30  # Solar mass to kg
        h = H0 / 100

        # Convert H(z) to SI units (1/s)
        H_z_si = H_z * 1000 / Mpc_to_m  # km/s/Mpc to 1/s

        # Critical density in SI units (kg/m^3)
        rho_crit_si = (3 * H_z_si**2) / (8 * np.pi * G)

        # Convert to 10^10 M_sun / (h^2 Mpc^3)
        rho_crit = rho_crit_si * Mpc_to_m**3 / Msun_to_kg * 1e-10 / h**2

        return rho_crit
    
    def H(self, z):
        """
        Get the Hubble parameter H at a given redshift.

        Parameters:
        -----------
        z : float
            Redshift.

        Returns:
        --------
        float
            The Hubble parameter H(z) in km/s/Mpc at redshift z.
        """
        return self._cosmo.hubble_parameter(z)
    
    def peak_height(self, M, z):
        """
        Calculate the peak-height delta_c / sigma(R(M), z).

        Parameters:
        -----------
        M : float
            Mass in solar masses over h.
        z : float
            Redshift.

        Returns:
        --------
        float
            The peak-height.
        """
        Om_z = self.Omega_m(z)
        delta_c_z = self.delta_c(Om_z)
        R = self.lagrangian_radius(M)
        sigma_R_z = self.sigma(R, z)
        return delta_c_z / sigma_R_z
    
    def lagrangian_radius(self, M):
        """
        Calculate the radius of the Lagrangian patch corresponding to a given mass M.

        Parameters:
        -----------
        M : float
            Mass in solar masses over h.

        Returns:
        --------
        float
            The Lagrangian radius in Mpc/h.
        """
        rho_m0 = self.critical_density(0) * self.Omega_m(0) * 1e10  # in M_sun / Mpc^3 * h^2
        R = (3 * M / (4 * np.pi * rho_m0))**(1/3)
        return R
    
    def sigma(self, R, z):
        """
        Calculate the RMS fluctuation sigma(R, z).

        Parameters:
        -----------
        R : float
            Radius in Mpc/h.
        z : float
            Redshift.

        Returns:
        --------
        float
            The RMS fluctuation sigma.
        """
        R = np.atleast_1d(R)
        sigma_squared = ssq_given_Dk(R[:, np.newaxis], self.k, self.Dk[z])
        return np.sqrt(sigma_squared)

    def dlnsigma_dlnR(self, R):
        """
        Calculate the derivative of the logarithm of sigma with respect to the logarithm of R.

        Parameters:
        -----------
        R : float
            Radius in Mpc/h.
        z : float
            Redshift.

        Returns:
        --------
        float
            The derivative dln(sigma)/dln(R).
        """
        R = np.atleast_1d(R)
        d_sigma_dR = d_s_given_Dk(R[:, np.newaxis], self.k, self.Dk[0])
        sigma_R = self.sigma(R, 0)

        # Avoiding newaxis to create a matrix when multiplying 
        return (R.flatten() / sigma_R) * d_sigma_dR
    
    def vfv(self, M, z, return_variables=False, halo_finder=best_fit_values_ROCKSTAR):
        """
        Calculate the multiplicity function, defined as νf(ν), where ν is the peak-height, as well as 
        related variables for a given mass and redshift.

        Parameters:
        -----------
        M : np.ndarray
            Mass of the halos in solar masses over h. 
        z : float
            Redshift. 
        return_variables : bool, optional
            If True, returns additional intermediate variables used in the calculation.
            If False, returns only the multiplicity function. Default is False.
        halo_finder : str, optional
            Descriptor for the best-fit values for the mass function parameters.
            Default is `ROCKSTAR`. Options are: `AHF`, `ROCKSTAR`, `SUBFIND`, and `VELOCIraptor`.

        Returns:
        --------
        np.ndarray or tuple
            If `return_variables` is False, returns an array of the multiplicity function νf(ν).
            If `return_variables` is True, returns a tuple containing:
            - M : Mass array
            - R : Lagrangian radius array corresponding to M
            - v : Peak-height array
            - dlnsdlnR : Derivative of log sigma with respect to log R
            - f(ν) : Multiplicity function values
        """
        if not isinstance(M, np.ndarray):
            raise RuntimeError("Masses should be an instance of numpy.ndarray")
        elif M.size < 2:
            raise RuntimeError("Masses size should be at least 2 (preferrebly much more).")
        R = self.lagrangian_radius(M)
        v = self.peak_height(M, z)
        dlnsdlnR = self.dlnsigma_dlnR(R)

        if return_variables:
            return M, R, v, dlnsdlnR, multiplicity_function(v, dlnsdlnR, self.Omega_m(z), best_fits[halo_finder])
        else:
            return multiplicity_function(v, dlnsdlnR, self.Omega_m(z), halo_finder)

    def dndlnM(self, M, z, halo_finder='ROCKSTAR'):
        """
        Calculate the halo mass function, dn/dlnM, which gives the number density of halos 
        per logarithmic mass interval.

        Parameters:
        -----------
        M : np.ndarray
            Mass of the halos in solar masses over h. 
        z : float
            Redshift. 
        halo_finder : str, optional
            Descriptor for the best-fit values for the mass function parameters.
            Default is `ROCKSTAR`. Options are: `AHF`, `ROCKSTAR`, `SUBFIND`, and `VELOCIraptor`.

        Returns:
        --------
        np.ndarray
            The halo mass function dn/dlnM, with units of number density per logarithmic mass interval.
        """
        
        M, R, v, dlnsdlnR, vfv = self.vfv(M, z, halo_finder=halo_finder, return_variables=True)

        return self.critical_density(0.) * self.Omega_m(0.) * 1e10 / M * vfv * np.gradient(np.log(v), np.log(M))
    

    def pbs_bias(self, M, z, halo_finder='ROCKSTAR', return_variables=False):
        """
        Calculate the PBS bias for halos of mass M at redshift z.

        Parameters:
        -----------
        M : np.ndarray
            Mass of the halos in solar masses over h. 
        z : float
            Redshift.
        halo_finder : str
            The halo finder used to determine the best-fit parameters. Default is 'ROCKSTAR'.
        return_variables : bool, optional
            If True, returns additional intermediate variables used in the calculation.
            If False, returns only the bias. Default is False.

        Returns:
        --------
        np.ndarray or tuple
            If `return_variables` is False, returns an array of the PBS bias.
            If `return_variables` is True, returns a tuple containing:
            - M : Mass array
            - R : Lagrangian radius array corresponding to M
            - v : Peak-height array
            - dlnsdlnR : Derivative of log sigma with respect to log R
            - vfv : Multiplicity function values
            - bias : PBS bias values
        """
        M, R, v, dlnsdlnR, vfv = self.vfv(M, z, halo_finder=halo_finder, return_variables=True)
        delta_c = self.delta_c(self.Omega_m(z))
        
        # Avoid division by zero or log of zero issues
        if np.any(v <= 0) or np.any(vfv <= 0):
            raise ValueError("Invalid values in v or vfv arrays; check input mass or redshift ranges.")
        
        bias = 1 - 1/delta_c * np.gradient(np.log(vfv), np.log(v))
        if return_variables:
            return M, R, v, dlnsdlnR, vfv, bias
        else:
            return bias

    def bias(self, M, z, halo_finder='ROCKSTAR'):
        """
        Calculate the corrected linear halo bias.

        Parameters:
        -----------
        M : np.ndarray
            Mass of the halos in solar masses over h. Can be a scalar or array.
        z : float
            Redshift.
        halo_finder : str
            The halo finder used to determine the best-fit parameters. Default is 'ROCKSTAR'.

        Returns:
        --------
        np.ndarray
            The corrected linear halo bias.
        """
        M, R, v, dlnsdlnR, vfv, b_pbs = self.pbs_bias(M, z, halo_finder=halo_finder, return_variables=True)
        Omz = self.Omega_m(z)
        S8 = self.sigma(8, 0.0) * np.sqrt(self.Omega_m(0.0)/0.3)
    
        return corrected_bias(b_pbs, Omz, dlnsdlnR, S8)

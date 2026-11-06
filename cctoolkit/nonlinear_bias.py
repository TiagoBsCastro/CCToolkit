"""
nonlinear_bias.py

Nonlinear halo-matter bias and halo auto power spectrum of DUCA III (Castro 2026).

The model predicts, for halos of a given mass at redshift z in a w0waCDM cosmology with
massive neutrinos, the cross bias b_hm(k) = P_hm(k) / P_mm(k) and the halo auto power
spectrum P_hh(k) (including the halo shot noise). It is calibrated on the DUCA
simulations for z < 2.5. All spectra refer to the cold matter (CDM plus baryons).

Two components are joined by a quintic smoothstep in ln k between 0.4 and 0.6 h/Mpc:

- the *base* model, calibrated at k <= 0.4 h/Mpc::

      b_hm(k) = b1 [1 + sum_mj w_mj phi_m B_j(x, y)],
      P_hh(k) = b_hm^2 P_mm [1 + sum_mj h_mj phi_m B_j(x, y)]
                + (1/n) [sum_m a_m phi_m + sum_mj a_mj phi_m C_j(x)],

  with x = k R_L, y = k / k_nl, and phi_m monomials of standardized halo, redshift and
  cosmology descriptors;
- the small-scale *extension*, which adds a compensated NFW one-halo term to the cross
  bias, uses bounded basis functions, and multiplies the departure from Poisson noise by
  an exclusion envelope. It is calibrated up to k = 1.13 (L), 1.61 (primary) and
  2.25 h/Mpc (high-resolution runs).

b1 = b_hm(k -> 0) is the Castro et al. (2025) peak-background-split bias times the
large-scale correction of Euclid Collaboration: Castro et al. (2024), recalibrated with
a cubic factor in the peak height.

The matter power spectrum P_mm entering the model is the nonlinear cold-matter spectrum;
the model was calibrated with the measured simulation spectrum. It is needed for P_hh
and, above k = 0.4 h/Mpc, for b_hm (the one-halo term is divided by it).

Domain of validity:

- z < 2.5, halos of at least 150 particles of the DUCA runs (for the Fiducial cosmology,
  M >~ 1.5e12 Msun/h from the high-resolution runs and M >~ 1.2e13 Msun/h from the
  primary ones), and the w0waCDM + massive-neutrino parameter ranges of the DUCA I
  simulations.
- The calibration uses narrow mass bins (Delta ln M = 0.2). In such bins the halo shot
  noise dominates P_hh at k >~ 0.3 h/Mpc (n b^2 P_mm <~ 0.2 there). For much denser
  samples, such as wide mass bins of low-mass halos, the small-scale clustering
  correction of P_hh is extrapolated: across the join the base model and its extension
  can then differ by tens of percent. b_hm is not affected.
- The join changes the form of the model across 0.4-0.6 h/Mpc; validate it before
  transforming the predictions to configuration space.

References
----------
- Castro (2026), "DUCA: Dynamic Universe Cosmological Analysis. III. Nonlinear halo bias
  and stochasticity in dynamical dark energy cosmologies".
- [arXiv:2504.07608](https://arxiv.org/abs/2504.07608): Castro et al. (2025), DUCA I.
- [arXiv:2409.01877](https://arxiv.org/abs/2409.01877): Euclid Collaboration: Castro
  et al. (2024), the linear halo bias.
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy.special import sici

from .utils import virial_Delta

PARAMETER_FILE = os.path.join(os.path.dirname(__file__), "data", "duca3_nonlinear_bias.json")

# Order of the descriptors theta in the parameter file.
DESCRIPTORS = (
    "peak_height",
    "z_scaled",
    "sigma_slope",
    "omega_cb_z",
    "S8",
    "w0",
    "wa",
    "mnu",
)

# Shape functions of x = k R_L and y = k / k_nl. Each vanishes as k^2 at k -> 0. The
# unbounded terms grow at most as x^2 or y^2; the bounded ones ("...b") saturate.
BASIS_TERMS = {
    "x2": lambda x, y: x**2,
    "x4": lambda x, y: x**4 / (1 + x**2),
    "x6": lambda x, y: x**6 / (1 + x**2) ** 2,
    "y2": lambda x, y: y**2,
    "y4": lambda x, y: y**4 / (1 + y**2),
    "y6": lambda x, y: y**6 / (1 + y**2) ** 2,
    "x2y2_x": lambda x, y: x**2 * y**2 / (1 + x**2),
    "x2b": lambda x, y: x**2 / (1 + x**2),
    "x4b": lambda x, y: x**4 / (1 + x**2) ** 2,
    "x6b": lambda x, y: x**6 / (1 + x**2) ** 3,
    "x8b": lambda x, y: x**8 / (1 + x**2) ** 4,
    "y2b": lambda x, y: y**2 / (1 + y**2),
    "y4b": lambda x, y: y**4 / (1 + y**2) ** 2,
    "y6b": lambda x, y: y**6 / (1 + y**2) ** 3,
    "y8b": lambda x, y: y**8 / (1 + y**2) ** 4,
    "x2y2b": lambda x, y: x**2 * y**2 / ((1 + x**2) * (1 + y**2)),
}

# Below these arguments the closed forms lose precision to cancellation; series are used.
_SERIES_TOP_HAT = 0.1
_SERIES_NFW = 1.0e-3
_SERIES_ENVELOPE = 1.0e-4


def smoothstep(k, lower, upper):
    """
    Quintic smoothstep in ln k: 0 at `lower` and below, 1 at `upper` and above.

    Parameters:
    -----------
    k : np.ndarray
        Wavenumbers in h/Mpc.
    lower, upper : float
        Wavenumbers where the step starts and ends, in h/Mpc.

    Returns:
    --------
    np.ndarray
        The weight of the small-scale extension.
    """
    t = np.clip(np.log(np.asarray(k, dtype=float) / lower) / np.log(upper / lower), 0.0, 1.0)
    return t**3 * (10.0 - 15.0 * t + 6.0 * t**2)


def top_hat_window(x):
    """Fourier transform of a spherical top hat, 3 (sin x - x cos x) / x^3."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < _SERIES_ENVELOPE
    x2 = x[small] ** 2
    out[small] = 1.0 - x2 / 10.0 + x2**2 / 280.0
    xb = x[~small]
    out[~small] = 3.0 * (np.sin(xb) - xb * np.cos(xb)) / xb**3
    return out


def one_minus_window(x, window="gaussian"):
    """1 - W(x) of the compensation window of the one-halo term, accurate as x -> 0."""
    x = np.abs(np.asarray(x, dtype=float))
    if window == "gaussian":
        return -np.expm1(-0.5 * x**2)
    if window == "lorentzian":
        return x**2 / (1.0 + x**2)
    if window != "top_hat":
        raise ValueError(f"Unknown compensation window '{window}'.")
    out = np.empty(x.shape)
    small = x < _SERIES_TOP_HAT
    x2 = x[small] ** 2
    out[small] = x2 / 10.0 - x2**2 / 280.0 + x2**3 / 15120.0
    xb = x[~small]
    out[~small] = 1.0 - 3.0 * (np.sin(xb) - xb * np.cos(xb)) / xb**3
    return out


def _nfw_mass(c):
    return np.log1p(c) - c / (1.0 + c)


def one_minus_nfw(k, virial_radius, concentration):
    """
    1 - u(k) of an NFW profile truncated at the virial radius, accurate as k -> 0.

    Parameters:
    -----------
    k : np.ndarray
        Wavenumbers in h/Mpc.
    virial_radius : np.ndarray
        Comoving virial radius in Mpc/h (broadcastable with k).
    concentration : np.ndarray
        Virial concentration (broadcastable with k).

    Returns:
    --------
    np.ndarray
        1 - u(k), with u the normalized Fourier transform of the profile.
    """
    k, rvir, c = np.broadcast_arrays(
        np.asarray(k, dtype=float),
        np.asarray(virial_radius, dtype=float),
        np.asarray(concentration, dtype=float),
    )
    x = k * rvir / c
    out = np.empty(x.shape)
    small = x < _SERIES_NFW
    if small.any():
        cs = c[small]
        t = 1.0 + cs
        moment = t**2 / 2.0 - 3.0 * t + 3.0 * np.log(t) + 1.0 / t + 1.5
        mean_square_radius = (rvir[small] / cs) ** 2 * moment / _nfw_mass(cs)
        out[small] = k[small] ** 2 * mean_square_radius / 6.0
    big = ~small
    xb, cb = x[big], c[big]
    si_outer, ci_outer = sici((1.0 + cb) * xb)
    si_inner, ci_inner = sici(xb)
    u = (
        np.sin(xb) * (si_outer - si_inner)
        - np.sin(cb * xb) / ((1.0 + cb) * xb)
        + np.cos(xb) * (ci_outer - ci_inner)
    ) / _nfw_mass(cb)
    out[big] = 1.0 - u
    return out


def concentration_duffy08(mass, redshift):
    """Duffy et al. (2008) virial concentration of the full sample (mass in Msun/h)."""
    return 7.85 * (np.asarray(mass, dtype=float) / 2.0e12) ** -0.081 * (
        1.0 + np.asarray(redshift, dtype=float)
    ) ** -0.71


class _Component:
    """One fitted smooth model (the base model or its small-scale extension)."""

    def __init__(self, payload, monomials, nu_centre):
        self.mean = np.asarray(payload["standardisation"]["mean"], dtype=float)
        self.std = np.asarray(payload["standardisation"]["std"], dtype=float)
        self.monomials = [tuple(m) for m in monomials]
        self.nu_centre = float(nu_centre)
        self.b1_coefficients = np.asarray(payload["b1_coefficients"], dtype=float)
        if self.b1_coefficients.size != 8:
            raise ValueError("Expected the cubic peak-height extension of b1 (8 coefficients).")
        self.cross_terms = list(payload["cross_basis"])
        self.auto_terms = list(payload["auto_basis"])
        self.stochastic_terms = list(payload["stochastic_basis"])
        nf = len(self.monomials)
        nc, na, ns = len(self.cross_terms), len(self.auto_terms), len(self.stochastic_terms)
        # Weights are stored feature-major: index m * n_terms + j.
        self.cross_weights = np.asarray(payload["cross_weights"], dtype=float).reshape(nf, nc)
        auto = np.asarray(payload["auto_weights"], dtype=float)
        if auto.size != nf * (na + 1 + ns):
            raise ValueError("Inconsistent number of auto-spectrum weights.")
        self.auto_weights = auto[: nf * na].reshape(nf, na)
        self.stochastic_amplitude = auto[nf * na : nf * (na + 1)]
        self.stochastic_weights = auto[nf * (na + 1) :].reshape(nf, ns)
        self.one_halo = payload.get("one_halo")
        self.exclusion = payload.get("exclusion")

    def features(self, theta):
        """Monomials of the standardized descriptors, shape (n_halo, n_features)."""
        z = (theta - self.mean) / self.std
        out = np.empty((theta.shape[0], len(self.monomials)))
        for column, monomial in enumerate(self.monomials):
            term = np.ones(theta.shape[0])
            for index in monomial:
                term = term * z[:, index]
            out[:, column] = term
        return out

    def b1(self, halos):
        """Large-scale bias: corrected PBS bias times a cubic factor in the peak height."""
        c = self.b1_coefficients
        s = halos["sigma_slope"]
        b1 = (
            halos["pbs_bias"]
            * c[0]
            * (1.0 + c[1] * halos["omega_cb_z"])
            * (1.0 + c[2] * s + c[3] * s**2)
            * (1.0 + c[4] * halos["S8"])
        )
        offset = halos["peak_height"] - self.nu_centre
        factor = 1.0 + c[5] * offset
        factor = factor + c[6] * offset**2
        factor = factor + c[7] * offset**3
        return b1 * factor

    @staticmethod
    def _sum(coefficients, terms, values):
        """sum_j coefficients[:, j] B_j, shape (n_halo, n_k); `values` maps names to B_j."""
        total = coefficients[:, 0, None] * values[terms[0]]
        for j in range(1, len(terms)):
            total += coefficients[:, j, None] * values[terms[j]]
        return total

    def evaluate(self, k, halos, theta, matter_power, want_auto):
        """
        Return (b_2h, b_hm, P_hh) on the (n_halo, n_k) grid; P_hh is None if not wanted.
        """
        phi = self.features(theta)
        radius = halos["lagrangian_radius"][:, None]
        x = k * radius
        y = k / halos["k_nl"][:, None]
        b1 = self.b1(halos)[:, None]
        names = self.cross_terms + (self.auto_terms + self.stochastic_terms if want_auto else [])
        values = {name: BASIS_TERMS[name](x, y) for name in dict.fromkeys(names)}
        b2h = b1 * (1.0 + self._sum(phi @ self.cross_weights, self.cross_terms, values))
        bhm = b2h
        if self.one_halo is not None:
            settings = self.one_halo
            omega = halos["omega_cb_z"][:, None]
            virial_radius = radius * np.cbrt(omega / virial_Delta(omega))
            concentration = settings["concentration_scale"] * concentration_duffy08(
                halos["mass"], halos["redshift"]
            )[:, None]
            bracket = one_minus_window(
                settings["compensation_scale"] * x, settings["compensation_window"]
            ) - one_minus_nfw(k, virial_radius, concentration)
            one_halo = settings["amplitude"] * 4.0 * np.pi * radius**3 / 3.0 * bracket
            bhm = b2h + one_halo / matter_power
        if not want_auto:
            return b2h, bhm, None
        clustering = b2h**2 * matter_power
        deterministic = self._sum(phi @ self.auto_weights, self.auto_terms, values)
        stochastic = (phi @ self.stochastic_amplitude)[:, None] + self._sum(
            phi @ self.stochastic_weights, self.stochastic_terms, values
        )
        inverse_density = 1.0 / halos["number_density"][:, None]
        if self.exclusion is None:
            halo_power = clustering + clustering * deterministic + inverse_density * stochastic
        else:
            q = k * self.exclusion["scale"] * radius
            envelope = top_hat_window(q) * np.exp(-0.5 * (self.exclusion["smoothing"] * q) ** 2)
            halo_power = (
                clustering
                + (1.0 - envelope) * inverse_density
                + clustering * deterministic
                + envelope * inverse_density * stochastic
            )
        return b2h, bhm, halo_power


class NonlinearHaloBias:
    """
    Evaluator of the DUCA III nonlinear halo bias and halo auto power spectrum.

    It works directly on the halo, redshift and cosmology descriptors. Most users will call
    `CosmologyCalculator.nonlinear_bias`, which computes the descriptors from the cosmology.

    Parameters:
    -----------
    parameter_file : str, optional
        JSON file with the calibrated parameters. Default: the file shipped with CCToolkit.

    Attributes:
    -----------
    blend : tuple
        Wavenumbers (h/Mpc) between which the base model is joined to its extension.
    k_max : dict
        Largest calibrated wavenumber (h/Mpc) by simulation set: the extension is
        calibrated up to 1.13 (L), 1.61 (primary) and 2.25 (high resolution).
    provenance : dict
        Hashes and commits of the files the parameters were exported from.
    """

    k_max = {"large_volume": 1.13, "primary": 1.61, "high_resolution": 2.25}
    z_max = 2.5

    def __init__(self, parameter_file=None):
        with open(parameter_file or PARAMETER_FILE) as handle:
            payload = json.load(handle)
        names = payload["descriptors"]
        expected = [
            "peak_height",
            "z_scaled",
            "logarithmic_sigma_slope",
            "omega_m",
            "s8",
            "w0",
            "wa",
            "neutrino_mass_ev",
        ]
        if names != expected:
            raise ValueError(f"Unexpected descriptors {names} in {parameter_file}.")
        self.blend = tuple(float(v) for v in payload["blend"])
        self.base = _Component(payload["base"], payload["features"], payload["nu_centre"])
        self.extension = _Component(
            payload["extension"], payload["features"], payload["nu_centre"]
        )
        self.provenance = payload.get("provenance", {})

    @staticmethod
    def _halos(n, **values):
        halos = {}
        for name, value in values.items():
            if value is None:
                halos[name] = None
                continue
            array = np.asarray(value, dtype=float)
            if array.ndim > 1:
                raise ValueError(f"'{name}' must be a scalar or a one-dimensional array.")
            halos[name] = np.broadcast_to(array, (n,)).astype(float)
        return halos

    def linear_bias(self, peak_height, sigma_slope, omega_cb_z, S8, pbs_bias):
        """
        The large-scale limit b1 = b_hm(k -> 0) of the model.

        Parameters:
        -----------
        peak_height : np.ndarray
            Peak height nu = delta_c / sigma(R_L, z) of the cold matter.
        sigma_slope : np.ndarray
            d ln sigma / d ln R at the Lagrangian radius.
        omega_cb_z : np.ndarray
            Cold-matter (CDM plus baryons) density parameter at the halo redshift.
        S8 : float or np.ndarray
            sigma8 of the cold matter times sqrt(Omega_m / 0.3), with the total present
            matter density Omega_m (including massive neutrinos).
        pbs_bias : np.ndarray
            Peak-background-split bias of the Castro et al. (2025) mass function.

        Returns:
        --------
        np.ndarray
            The linear bias b1.
        """
        n = np.broadcast(peak_height, sigma_slope, omega_cb_z, S8, pbs_bias).size
        halos = self._halos(
            n,
            peak_height=peak_height,
            sigma_slope=sigma_slope,
            omega_cb_z=omega_cb_z,
            S8=S8,
            pbs_bias=pbs_bias,
        )
        return self.base.b1(halos)

    def predict(
        self,
        k,
        lagrangian_radius,
        k_nl,
        peak_height,
        redshift,
        sigma_slope,
        omega_cb_z,
        S8,
        w0,
        wa,
        mnu,
        pbs_bias,
        mass=None,
        number_density=None,
        matter_power=None,
    ):
        """
        Cross bias b_hm(k) and halo auto power spectrum P_hh(k).

        Every halo argument is a scalar or an array of length n_halo (one entry per halo
        mass or mass bin); scalars are broadcast.

        Parameters:
        -----------
        k : np.ndarray
            Wavenumbers in h/Mpc: one-dimensional (shared by all halos) or (n_halo, n_k).
        lagrangian_radius : np.ndarray
            R_L = (3 M / 4 pi rho_cb)^(1/3) in Mpc/h, with rho_cb the mean cold-matter density.
        k_nl : np.ndarray
            Wavenumber (h/Mpc) where the linear cold-matter Delta^2(k) = 1 at the redshift.
        peak_height, sigma_slope, omega_cb_z, S8, pbs_bias : np.ndarray
            As in `linear_bias`.
        redshift : np.ndarray
            Redshift (calibrated for z < 2.5).
        w0, wa : np.ndarray
            CPL dark-energy parameters.
        mnu : np.ndarray
            Sum of the neutrino masses in eV.
        mass : np.ndarray, optional
            Halo (mean) mass in Msun/h. Needed above k = 0.4 h/Mpc (one-halo term).
        number_density : np.ndarray, optional
            Number density of the halo sample in (h/Mpc)^3. Needed for P_hh.
        matter_power : np.ndarray, optional
            Nonlinear cold-matter power spectrum in (Mpc/h)^3, broadcastable to the
            (n_halo, n_k) grid. Needed for P_hh, and for b_hm above k = 0.4 h/Mpc.

        Returns:
        --------
        tuple
            (b_hm, P_hh), each of shape (n_halo, n_k); P_hh is None unless both
            `number_density` and `matter_power` are given.
        """
        n = np.broadcast(
            lagrangian_radius, k_nl, peak_height, redshift, sigma_slope, omega_cb_z, S8,
            w0, wa, mnu, pbs_bias,
        ).size
        halos = self._halos(
            n,
            lagrangian_radius=lagrangian_radius,
            k_nl=k_nl,
            peak_height=peak_height,
            redshift=redshift,
            sigma_slope=sigma_slope,
            omega_cb_z=omega_cb_z,
            S8=S8,
            w0=w0,
            wa=wa,
            mnu=mnu,
            pbs_bias=pbs_bias,
            mass=mass,
            number_density=number_density,
        )
        k = np.asarray(k, dtype=float)
        if np.any(k <= 0):
            raise ValueError("Wavenumbers must be positive.")
        if k.ndim == 1:
            # Shared wavenumbers: the join weight depends on k only.
            weight = smoothstep(k, *self.blend)[None, :]
            k = np.broadcast_to(k, (n, k.size))
        elif k.ndim == 2 and k.shape[0] == n:
            weight = smoothstep(k, *self.blend)
        else:
            raise ValueError("k must be one-dimensional or of shape (n_halo, n_k).")
        want_auto = number_density is not None and matter_power is not None
        if matter_power is not None:
            matter_power = np.broadcast_to(np.asarray(matter_power, dtype=float), k.shape)
        if np.any(weight > 0):
            if matter_power is None:
                raise ValueError(
                    f"Above k = {self.blend[0]} h/Mpc the cross bias needs the matter power spectrum."
                )
            if halos["mass"] is None:
                raise ValueError(f"Above k = {self.blend[0]} h/Mpc the halo mass is needed.")
        theta = np.column_stack(
            [
                halos["peak_height"],
                halos["redshift"] / (1.0 + halos["redshift"]),
                halos["sigma_slope"],
                halos["omega_cb_z"],
                halos["S8"],
                halos["w0"],
                halos["wa"],
                halos["mnu"],
            ]
        )
        # Evaluate each component only on the wavenumber columns where it contributes, and
        # blend the two on the columns of the join.
        base_cols = np.flatnonzero(np.any(weight < 1.0, axis=0))
        ext_cols = np.flatnonzero(np.any(weight > 0.0, axis=0))
        bias = np.empty(k.shape)
        halo = np.empty(k.shape) if want_auto else None
        results = []
        for component, cols in ((self.base, base_cols), (self.extension, ext_cols)):
            if cols.size == 0:
                results.append((None, None))
                continue
            power = None if matter_power is None else matter_power[:, cols]
            _, bhm, phh = component.evaluate(k[:, cols], halos, theta, power, want_auto)
            bias[:, cols] = bhm
            if want_auto:
                halo[:, cols] = phh
            results.append((bhm, phh))
        join, in_base, in_ext = np.intersect1d(base_cols, ext_cols, return_indices=True)
        if join.size:
            w = weight[:, join]
            (base_b, base_p), (ext_b, ext_p) = results
            bias[:, join] = (1.0 - w) * base_b[:, in_base] + w * ext_b[:, in_ext]
            if want_auto:
                halo[:, join] = (1.0 - w) * base_p[:, in_base] + w * ext_p[:, in_ext]
        return bias, halo

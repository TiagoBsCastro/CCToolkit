import json
import os

import numpy as np
import pytest

from cctoolkit import CosmologyCalculator, NonlinearHaloBias
from cctoolkit.nonlinear_bias import smoothstep

REFERENCE = os.path.join(os.path.dirname(__file__), "data", "duca3_reference.json")

PARAMS = {
    "H0": 67.321,
    "Ob0": 0.0494,
    "Om0": 0.3158,
    "sigma8": 0.8102,
    "ns": 0.9661,
    "TCMB": 2.7255,
    "mnu": 0.06,
    "num_massive_neutrinos": 1,
    "w0": -0.9,
    "wa": -0.3,
}


@pytest.fixture(scope="module")
def model():
    return NonlinearHaloBias()


@pytest.fixture(scope="module")
def reference():
    with open(REFERENCE) as handle:
        rows = json.load(handle)
    rows.pop("repository_commit")
    return {name: np.asarray(values) for name, values in rows.items()}


@pytest.fixture(scope="module")
def calculator():
    return CosmologyCalculator(PARAMS, zmax=2.5)


def _descriptor_args(rows):
    return (
        rows["lagrangian_radius"],
        rows["k_nl"],
        rows["peak_height"],
        rows["redshift"],
        rows["logarithmic_sigma_slope"],
        rows["omega_m"],
        rows["s8"],
        rows["w0"],
        rows["wa"],
        rows["neutrino_mass_ev"],
        rows["pbs_bias"],
    )


# The reference rows were evaluated with the code that calibrated the model
# (github.com/TiagoBsCastro/NonLinearBias), across the calibrated range and the join.
def test_reproduces_reference_rows(model, reference):
    bias, halo = model.predict(
        reference["k"][:, None],
        *_descriptor_args(reference),
        mass=reference["mean_mass"],
        number_density=reference["number_density"],
        matter_power=reference["matter_power"][:, None],
    )
    np.testing.assert_allclose(bias[:, 0], reference["b_hm"], rtol=1e-12)
    np.testing.assert_allclose(halo[:, 0], reference["P_hh"], rtol=1e-12)


def test_shared_wavenumbers_match_row_by_row(model, reference):
    rows = {name: values[:50] for name, values in reference.items()}
    power = 1.0e3 * (rows["k"] / 0.1) ** -1.4 + 2.0  # the reference P_mm(k)
    bias, halo = model.predict(
        rows["k"],
        *_descriptor_args(rows),
        mass=rows["mean_mass"],
        number_density=rows["number_density"],
        matter_power=power,
    )
    assert bias.shape == halo.shape == (50, 50)
    np.testing.assert_allclose(np.diag(bias), rows["b_hm"], rtol=1e-12)
    np.testing.assert_allclose(np.diag(halo), rows["P_hh"], rtol=1e-12)


def test_large_scale_limit(model):
    peak = np.array([0.8, 2.0, 3.5])
    args = dict(sigma_slope=-0.6, omega_cb_z=0.5, S8=0.85, pbs_bias=np.array([0.9, 2.0, 6.0]))
    b1 = model.linear_bias(peak, **args)
    bias, _ = model.predict(
        np.array([1e-4]), np.array([1.5, 5.0, 12.0]), 0.3, peak, 0.5,
        args["sigma_slope"], args["omega_cb_z"], args["S8"], -1.0, 0.0, 0.06, args["pbs_bias"],
    )
    np.testing.assert_allclose(bias[:, 0], b1, rtol=1e-5)


def test_join_is_continuous(model):
    k = np.geomspace(0.3, 0.8, 2001)
    power = 1.0e3 * (k / 0.1) ** -1.4 + 2.0
    bias, halo = model.predict(
        k, 4.0, 0.3, 2.2, 0.5, -0.6, 0.5, 0.85, -1.0, 0.0, 0.06, 2.5,
        mass=1e14, number_density=1e-5, matter_power=power,
    )
    for values in (bias[0], halo[0]):
        steps = np.abs(np.diff(np.log(values)))
        assert steps.max() < 1e-3
    assert smoothstep(0.4, 0.4, 0.6) == 0.0 and smoothstep(0.6, 0.4, 0.6) == 1.0


def test_small_scales_need_matter_power(model):
    with pytest.raises(ValueError):
        model.predict(np.array([0.5]), 4.0, 0.3, 2.2, 0.5, -0.6, 0.5, 0.85, -1.0, 0.0, 0.06, 2.5,
                      mass=1e14)


def test_calculator_descriptors(calculator):
    masses = np.geomspace(1e12, 1e15, 7)
    d = calculator.halo_descriptors(masses, 0.5)
    np.testing.assert_allclose(d["lagrangian_radius"], calculator.lagrangian_radius(masses))
    np.testing.assert_allclose(d["peak_height"], calculator.peak_height(masses, 0.5))
    assert np.all(np.diff(d["pbs_bias"]) > 0)
    assert 0.0 < d["k_nl"] < 5.0 and 0.5 < d["S8"] < 1.2
    assert d["omega_cb_z"] == pytest.approx(calculator.Omega_m(0.5))
    # The bin-averaged bias lies between the point values at the bin edges.
    edges = np.geomspace(1e13, 1e14, 3)
    binned = calculator.halo_descriptors(edges, 0.5, bins=True)
    points = calculator.halo_descriptors(edges, 0.5)
    assert np.all(binned["pbs_bias"] > points["pbs_bias"][:-1])
    assert np.all(binned["pbs_bias"] < points["pbs_bias"][1:])
    assert np.all((binned["mass"] > edges[:-1]) & (binned["mass"] < edges[1:]))
    assert np.all(binned["number_density"] > 0)


def test_calculator_nonlinear_bias(calculator):
    k = np.geomspace(1e-3, 1.5, 200)
    power = 2.0e4 * (k / 0.02) / (1.0 + (k / 0.02) ** 2.4) + 50.0
    masses = np.geomspace(1e12, 1e15, 4)
    bias, halo = calculator.nonlinear_bias(k[k < 0.4], masses, 0.5)
    assert bias.shape == (4, np.sum(k < 0.4)) and halo is None
    bias, halo, d = calculator.nonlinear_bias(
        k, np.geomspace(1e12, 1e15, 5), 1.0, matter_power=power, bins=True, return_variables=True
    )
    assert bias.shape == halo.shape == (4, k.size)
    assert np.all(np.isfinite(bias)) and np.all(halo > 0)
    b1 = NonlinearHaloBias().linear_bias(
        d["peak_height"], d["sigma_slope"], d["omega_cb_z"], d["S8"], d["pbs_bias"]
    )
    np.testing.assert_allclose(bias[:, 0], b1, rtol=1e-3)
    # On large scales the auto spectrum is clustering plus stochastic noise near 1/n.
    noise = halo[:, 0] - bias[:, 0] ** 2 * power[0]
    assert np.all(np.abs(noise * d["number_density"] - 1.0) < 0.5)

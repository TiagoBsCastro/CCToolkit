import matplotlib.pyplot as plt
import numpy as np
from cctoolkit.cosmology import CosmologyCalculator
try:
    from ConceptSpectra import get_power, get_modes
    skip = False
except:
    skip = True

if not skip:

    # Define the cosmological parameters for a model with a dark energy fluid for concept spectra
    params = {
        # Cosmology parameters (with no cosmological constant; DE fluid is used instead)
        "H0": 67,                  # Hubble constant in km/s/Mpc
        "Omega_b": 0.049,          # Baryon density parameter
        "Omega_cdm": 0.27,         # Cold dark matter density parameter
        "Omega_Lambda": 0,         # No cosmological constant since DE fluid is used
        "w0_fld": -0.9,            # Dark energy fluid equation-of-state parameter, w0
        "wa_fld": 0.2,             # Dark energy fluid evolution parameter, wa
        "cs2_fld": 1e-7,           # Sound speed squared for the dark energy fluid
        "use_ppf": "no",           # Use ppf or fluid DE description
        # Primordial spectrum parameters
        "A_s": 2.1e-9,             # Amplitude of primordial perturbations
        "n_s": 0.96,               # Spectral index
        # Precision parameters for CLASS
        "l_max_g": 100,
        "l_max_pol_g": 100,
        "radiation_streaming_approximation": 3,
        "l_max_ur": 100,
        "ur_fluid_approximation": 3,
        "evolver": 0,
        "recfast_Nz0": 1e5,
        "tol_thermo_integration": 1e-6,
        "perturb_sampling_stepsize": 0.01,
        # Output settings
        "output": "dTk",
        # Generate modes as a formatted string using CO*N*CEPTSpectra's get_modes
        "k_output_values": get_modes(1e-3, 1e1, 30, as_str=True)
    }

    # Set the scale factor at which to compute the power spectrum.
    # Here we choose a = 1.0, but this can be modified as needed.
    a   = 1.0
    z   = 1/a - 1
    zta = (1+z)/(1/2)**(2/3) - 1
    ata = 1/(1+zta)

    # Compute the DE fluid power spectrum using the "fld" species.
    # Using the "nbody" gauge transformation for an example.
    modes, de_power = get_power(params, "fld", "nbody", a=ata)

    # Params for cctoolkit
    params = {
            "H0": params['H0'],
            "Ob0": params['Omega_b'],
            "Om0": params['Omega_cdm'] + params['Omega_b'],
            "As": params['A_s'],
            "ns": params['n_s'],
            "TCMB": 2.7255,
            "mnu": 0.00,
            "num_massive_neutrinos": 0,
            "w0": params["w0_fld"],
            "wa": params["wa_fld"],
        }
    cc = CosmologyCalculator(params)
    masses = np.logspace(13, 15.5, num=100)
    hmf_castro25  = cc.dndlnM(masses, 0, model="castro25")
    hmf_castro25b = cc.dndlnM(masses, 0, model="castro25", PkDE=np.transpose([modes, de_power]))

    plt.plot(masses, hmf_castro25b/hmf_castro25-1)

    data = np.loadtxt("data_plot_model_w_0=-0.9_w_a=0.2_cs2=1.0e-07_z=0.0.txt")
    plt.plot(data[:, 0], data[:, 1])
    plt.xscale('log')
    plt.ylim([-.1,.1])
    plt.savefig("cde.pdf")
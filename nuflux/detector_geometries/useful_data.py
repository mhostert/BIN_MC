import numpy as np

from scipy import interpolate

from helpers import material, subs, comp, unif


def cs_interp():
    """Returns necessary cross neutrino cross sections from xsecs"""
    log10E, sigmaeo, sigmamuo, _, sigmaebaro, sigmamubaro, _ = np.genfromtxt(
        "xsecs/XCC.dat", unpack=True
    )
    exs = 10 ** (log10E)
    exs = np.append(exs, 1e4)

    # adding values linearily up to 10 TeV
    sigmae = np.append(sigmaeo, 7.046434082801659e-1)
    sigmamu = np.append(sigmamuo, 7.046434082801659e-1)
    sigmaebar = np.append(sigmaebaro, 3.758631195493489e-1)
    sigmamubar = np.append(sigmamubaro, 3.758631195493489e-1)

    # generating interpolators
    sigmanue = interpolate.interp1d(
        exs, sigmae * exs * 1e-38, bounds_error=False, fill_value=0.0
    )
    sigmanuebar = interpolate.interp1d(
        exs, sigmaebar * exs * 1e-38, bounds_error=False, fill_value=0.0
    )
    sigmanumu = interpolate.interp1d(
        exs, sigmamu * exs * 1e-38, bounds_error=False, fill_value=0.0
    )
    sigmanumubar = interpolate.interp1d(
        exs, sigmamubar * exs * 1e-38, bounds_error=False, fill_value=0.0
    )

    return sigmanue, sigmanuebar, sigmanumu, sigmanumubar


# density in g/cm**3; atomic mass in g/mol
Si = subs(2.329, 28.0855, 28)
WSi2 = subs(9.3, 240.01, 240)
Fe = subs(7.874, 55.845, 56)
Al = subs(2.7, 26.981539, 27)
W = subs(19.3, 183.84, 184)
Cu = subs(8.96, 63.546, 64)
PS = subs(1.05, 104.1, 104)

# from CLICdet paper
hcal_CLICdet = comp(
    [[Fe, 20 / 26.5], [Al, 0.7 / 26.5], [Cu, 0.1 / 26.5], [PS, 3 / 26.5]]
)
ecal_CLICdet = comp([[W, 1.9 / 5.05], [Cu, 2.3 / 5.05], [Si, 0.5 / 5.05]])

# from online
EARTH_DENSITY = unif(5.51).N

# Parameter sets; from muTristan (Hamada et al., 2022), and MuCoL (Intl Muon Coll Collab, 2023)
COMMON_dTHETA = 5.88e-4  # radians

mutristan_small = {
    "name": "μTRISTAN (s)",
    "beam_p0": 1e3,
    "pmax": 3e3,
    "pmin": 0,
    "C": 3e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 1.71e-4,
    "circular": False,
    "Nmu": (1.4e10 * 3.6/2.3 * (1-1/np.e) * 40) * (50) * (365.25 * 24 * 3600 * 0.7), #3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "syr": 365.25 * 24 * 3600 * 0.7,
    "bunch": 40,
    "finj": 50,
}
mutristan_large = {
    "name": "μTRISTAN (L)",
    "beam_p0": 3e3,
    "pmax": 9e3,
    "pmin": 0,
    "C": 9e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 99e-5,
    "circular": False,
    "Nmu": 3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "syr": 365.25 * 24 * 3600 * 0.7,
    "bunch": 40,
    "finj": 50,
}
mucol_s1 = {
    "name": "MuCoL (s1)",
    "beam_p0": 1.5e3,
    "pmax": 4.5e3,
    "pmin": 0,
    "C": 4.5e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
    "circular": False,
    "Nmu": 4.9e9 * 4.5e3 * 1.2e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 5,
}
mucol_s2 = {
    "name": "IMCC-II",
    "beam_p0": 5e3,
    "pmax": 15e3,
    "pmin": 0,
    "C": 10e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
    "circular": False,
    "Nmu": (1.8e12 * 2) * (5) * (1.2e7), #1.8e9 * 10e3 * 1.2e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 5,
}
scd_cern = {
    "name": "SCD CERN",
    "beam_p0": 5e3,
    "pmax": 15e3,
    "pmin": 0,
    "C": 866700,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 5.90e-4,
    "circular": False,
    "Nmu": 1.8e9 * 10e3 * 139 / 365 * 3.154e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 5,
}
mokhov = {
    "name": "Mohkov et al. (Fermilab)",
    "beam_p0": 75e1,
    "pmax": 75e1 * 3,
    "pmin": 0,
    "C": 273000,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 1.05e-3,
    "circular": False,
    "Nmu": (2e12 * 2) * (15) * (200 * 24 *3600),# 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 15,
}

parameters = {
    "mutristan s": mutristan_small,
    "mutristan l": mutristan_large,
    "mucol s1": mucol_s1,
    "mucol s2": mucol_s2,
    "scd cern": scd_cern,
    "mokhov": mokhov,
}


# These are weights to be added to each event from GENIE based on the particle that generated it and the detector in which it interacted. Necessary since I've generated the same amount of detector interactions in each region but they don't all have the same weight. These weights are the percentages of the total count these components have based on large simulations.
wmuTs = {
    "Name": "μTRISTAN (s)",
    "tc": 8.32e10,
    "nue": {"MD": 52.8, "SB": 1.6, "SM": 2.3, "HC": 6.2, "EC": 0.0, "NO": 0.1},
    "numubar": {"MD": 30.9, "SB": 0.9, "SM": 1.4, "HC": 3.6, "EC": 0.0, "NO": 0.1},
}
wmokhov = {
    "Name": "Mohkov et al. (Fermilab)",
    "tc": 8.94e10,
    "nue": {"MD": 15.1, "SB": 0.3, "SM": 0.5, "HC": 7.5, "EC": 1.4, "NO": 5.8},
    "numubar": {"MD": 8.9, "SB": 0.2, "SM": 0.3, "HC": 4.4, "EC": 0.8, "NO": 3.5},
    "nuebar": {"MD": 7.7, "SB": 0.2, "SM": 0.3, "HC": 3.8, "EC": 0.7, "NO": 3.0},
    "numu": {"MD": 17.6, "SB": 0.4, "SM": 0.6, "HC": 8.7, "EC": 1.6, "NO": 6.8},
}
wmucols2 = {
    "Name": "IMCC-II",
    "tc": 1.70e11,
    "nue": {"MD": 14.5, "SB": 0.3, "SM": 0.5, "HC": 7.1, "EC": 0.9, "NO": 7.0},
    "nuebar": {"MD": 7.5, "SB": 0.2, "SM": 0.3, "HC": 3.7, "EC": 0.5, "NO": 3.6},
    "numu": {"MD": 16.9, "SB": 0.4, "SM": 0.6, "HC": 8.3, "EC": 1.0, "NO": 8.3},
    "numubar": {"MD": 8.8, "SB": 0.2, "SM": 0.3, "HC": 4.3, "EC": 0.5, "NO": 4.3},
}

weights = {"muTs": wmuTs, "mokhov": wmokhov, "mucols2": wmucols2}

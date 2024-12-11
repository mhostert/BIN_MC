import numpy as np

colls_types_to_beam_cases = {
    "mu+e-": [["nue", "left"], ["numubar", "left"]],
    "mu+mu+": [
        ["nue", "left"],
        ["numubar", "left"],
        ["nue", "right"],
        ["numubar", "right"],
    ],
    "mu+mu-": [
        ["nue", "left"],
        ["numubar", "left"],
        ["nuebar", "right"],
        ["numu", "right"],
    ],
}
acc_colls_types = ["mu+e-", "mu+mu+", "mu+mu-"]
acc_colls_dict = {
    "mu+e-": r"$\mu^+ e^-$",
    "mu+mu+": r"$\mu^+ \mu^+$",
    "mu+mu-": r"$\mu^+ \mu^-$",
}

# Parameter sets; from muTristan (Hamada et al., 2022), and MuCoL (Intl Muon Coll Collab, 2023)
COMMON_dTHETA = 5.88e-4  # radians

##########################################
mut_2tev = {
    "name": r"$\mu$TRISTAN 2 TeV",
    "short_name": "mut_2tev",
    "beam_p0": 1e3,
    "C": 3e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": 1.70e-4,
<<<<<<< HEAD
    "Nmu_per_bunch": (
        1.4e10 * 3.6 / 2.3 * (1 - 1 / np.e)
    ),  # 3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "duty_factor": 1 / np.pi,
=======
    "circular": False,
    "Nmu_per_bunch": (
        1.4e10 * 3.6 / 2.3 * (1 - 1 / np.e)
    ),  # 3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "duty_factor": 0.7,
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
    "bunch_multiplicity": 40,
    "finj": 50,
    "Lss": 75e2,
    "collision_type": "mu+mu+",
    "muon_polarization": 0.0,
}

##########################################
mut_2tev_pos_pol = mut_2tev.copy()
mut_2tev_pos_pol["name"] = r"$\mu$TRISTAN 2 TeV ($P=0.8$)"
mut_2tev_pos_pol["short_name"] = "mut_2tev_pol"
mut_2tev_pos_pol["muon_polarization"] = 0.8

##########################################
mut_2tev_neg_pol = mut_2tev.copy()
mut_2tev_neg_pol["name"] = r"$\mu$TRISTAN 2 TeV ($P=-0.8$)"
mut_2tev_neg_pol["short_name"] = "mut_2tev_pol"
mut_2tev_neg_pol["muon_polarization"] = -0.8

##########################################
mut_6tev = {
    "name": r"$\mu$TRISTAN 6 TeV",
    "short_name": "mut_6tev",
    "beam_p0": 3e3,
    "C": 9e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 99e-5,
<<<<<<< HEAD
    "Nmu_per_bunch": 3.6e-9 / (1.6e-19) * (1 - 1 / np.e),
    "duty_factor": 1 / np.pi,
=======
    "circular": False,
    "Nmu_per_bunch": 3.6e-9 / (1.6e-19) * (1 - 1 / np.e),
    "duty_factor": 0.7,
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
    "bunch_multiplicity": 40,
    "finj": 50,
    "Lss": 75e2,
    "collision_type": "mu+mu+",
    "muon_polarization": 0.0,
}

##########################################
muc_3tev = {
    "name": r"MuC 3 TeV",
    "short_name": "muc_3tev",
    "beam_p0": 1.5e3,
    "C": 4.5e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
<<<<<<< HEAD
    "Nmu_per_bunch": 1.8e12,  # 4.9e9 * 4.5e3 * 1.2e7,
    "duty_factor": 1 / np.pi,
=======
    "circular": False,
    "Nmu_per_bunch": (1.8e12 * 2),  # 4.9e9 * 4.5e3 * 1.2e7,
    "duty_factor": 1.2 / np.pi,
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
    "bunch_multiplicity": 1,
    "finj": 5,
    "Lss": 100e2,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.0,
}

##########################################
muc_10tev = {
    "name": r"MuC 10 TeV",
    "short_name": "muc_10tev",
    "beam_p0": 5e3,
    "C": 10e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
<<<<<<< HEAD
    "Nmu_per_bunch": 1.8e12,  # 1.8e9 * 10e3 * 1.2e7, #NOTE: Why factor of 2 here?
    "duty_factor": 1 / np.pi,
=======
    "circular": False,
    "Nmu_per_bunch": (1.8e12 * 2),  # 1.8e9 * 10e3 * 1.2e7, #NOTE: Why factor of 2 here?
    "duty_factor": 1.2 / np.pi,
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
    "bunch_multiplicity": 1,
    "finj": 5,
    "Lss": 100e2,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.0,
}

##########################################
muc_30tev = {
    "name": r"MuC 30 TeV",
    "short_name": "muc_30tev",
    "beam_p0": 15e3,
    "C": 10e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
<<<<<<< HEAD
    "Nmu_per_bunch": 1.8e12,  # 1.8e9 * 10e3 * 1.2e7,
    "duty_factor": 1 / np.pi,
=======
    "circular": False,
    "Nmu_per_bunch": (1.8e12 * 2),  # 1.8e9 * 10e3 * 1.2e7,
    "duty_factor": 1.2 / np.pi,
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
    "bunch_multiplicity": 1,
    "finj": 5,
    "Lss": 100e2,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.0,
}

##########################################
muc_1p5tev = {
    "name": r"MuC 1.5 TeV",
    "short_name": "muc_1p5tev",
    "beam_p0": 750,
    "C": 273000,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 5.93e-4,
<<<<<<< HEAD
    "Nmu_per_bunch": 1.8e12,  # 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "duty_factor": 1 / np.pi,
=======
    "circular": False,
    "Nmu_per_bunch": (1.8e12 * 2),  # 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "duty_factor": 1.2 / np.pi,
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
    "bunch_multiplicity": 1,
    "finj": 15,
    "Lss": 50e2,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.0,
<<<<<<< HEAD
}

##########################################
nf_100gev_pos = {
    "name": r"NF 50 GeV ($P=1$)",
    "short_name": "muc_1p5tev",
    "beam_p0": 100,
    "C": 1000e2,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 5.93e-4,
    "Nmu_per_bunch": 1.8e12,  # 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "duty_factor": 1 / np.pi,
    "bunch_multiplicity": 1,
    "finj": 15,
    "Lss": 100e2,
    "collision_type": "mu+mu-",
    "muon_polarization": 1.0,
}
nf_100gev_neg = {
    "name": r"NF 50 GeV ($P=1$)",
    "short_name": "muc_1p5tev",
    "beam_p0": 100,
    "C": 1000e2,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 5.93e-4,
    "Nmu_per_bunch": 1.8e12,  # 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "duty_factor": 1 / np.pi,
    "bunch_multiplicity": 1,
    "finj": 15,
    "Lss": 100e2,
    "collision_type": "mu+mu-",
    "muon_polarization": -1.0,
=======
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
}

##########################################
# scd_cern = {
#     "name": r"SCD CERN",
#     "beam_p0": 5e3,
#     "C": 866700,
#     "beam_dpop": 1e-3,
#     "beam_dtheta": COMMON_dTHETA,  # 5.88e-4,
#     "Nmu": 1.8e9 * 10e3 * 139 / 365 * 3.154e7,
#     "syr": 1.2e7,
#     "bunch": 1,
#     "finj": 5,
#     "Lss": 50,
#     "collision_type": "mu+mu-",
#     "muon_polarization": 0.0,
# }

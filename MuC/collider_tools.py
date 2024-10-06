import numpy as np

colls_types_to_part = {
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
    "circular": False,
    "Nmu": (1.4e10 * 3.6 / 2.3 * (1 - 1 / np.e) * 40)
    * (50)
    * (
        365.25 * 24 * 3600 * 0.7
    ),  # 3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "syr": 365.25 * 24 * 3600 * 0.7,
    "bunch": 40,
    "finj": 50,
    "Lss": 75,
    "collision_type": "mu+mu+",
    "muon_polarization": 0.5,
}

##########################################
mut_2tev_pol = mut_2tev.copy()
mut_2tev_pol["name"] = "$\mu$TRISTAN 2 TeV (P=0.8)"
mut_2tev_pol["short_name"] = "mut_2tev_pol"
mut_2tev_pol["muon_polarization"] = 0.8

##########################################
mut_6tev = {
    "name": r"$\mu$TRISTAN 6 TeV",
    "short_name": "mut_6tev",
    "beam_p0": 3e3,
    "C": 9e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 99e-5,
    "circular": False,
    "Nmu": 3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "syr": 365.25 * 24 * 3600 * 0.7,
    "bunch": 40,
    "finj": 50,
    "Lss": 75,
    "collision_type": "mu+mu+",
    "muon_polarization": 0.5,
}

##########################################
muc_3tev = {
    "name": r"MuC 3 TeV",
    "short_name": "muc_3tev",
    "beam_p0": 1.5e3,
    "C": 4.5e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
    "circular": False,
    "Nmu": 4.9e9 * 4.5e3 * 1.2e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 5,
    "Lss": 200,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.5,
}

##########################################
muc_10tev = {
    "name": r"MuC 10 TeV",
    "short_name": "muc_10tev",
    "beam_p0": 5e3,
    "C": 10e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
    "circular": False,
    "Nmu": (1.8e12 * 2) * (5) * (1.2e7),  # 1.8e9 * 10e3 * 1.2e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 5,
    "Lss": 200,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.5,
}

##########################################
muc_30tev = {
    "name": r"MuC 30 TeV",
    "short_name": "muc_30tev",
    "beam_p0": 15e3,
    "C": 10e5,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,
    "circular": False,
    "Nmu": (1.8e12 * 2) * (5) * (1.2e7),  # 1.8e9 * 10e3 * 1.2e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 5,
    "Lss": 200,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.5,
}

##########################################
muc_1p5tev = {
    "name": r"MuC 1.5 TeV",
    "short_name": "muc_1p5tev",
    "beam_p0": 750,
    "C": 273000,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 5.93e-4,
    "circular": False,
    "Nmu": (2e12 * 2)
    * (15)
    * (200 * 24 * 3600),  # 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 15,
    "Lss": 50,
    "collision_type": "mu+mu-",
    "muon_polarization": 0.5,
}

##########################################
# scd_cern = {
#     "name": r"SCD CERN",
#     "beam_p0": 5e3,
#     "C": 866700,
#     "beam_dpop": 1e-3,
#     "beam_dtheta": COMMON_dTHETA,  # 5.88e-4,
#     "circular": False,
#     "Nmu": 1.8e9 * 10e3 * 139 / 365 * 3.154e7,
#     "syr": 1.2e7,
#     "bunch": 1,
#     "finj": 5,
#     "Lss": 50,
#     "collision_type": "mu+mu-",
#     "muon_polarization": 0.5,
# }

##############################################################################################################################
# NOTE: I am translating all this by hand info into functions.... give me a few days
##############################################################################################################################

"""
 These are weights to be added to each event from GENIE based on the particle that generated it and the detector in which it interacted.
 Necessary since I've generated the same amount of detector interactions in each region but they don't all have the same weight.
 These weights are the percentages of the total count these components have based on large simulations.
 """
wmuTs = {
    "Name": "$\mu$TRISTAN (s)",
    "tc": 81264561993.55392,
    "nue": {
        "MD": 52.75018261623187,
        "HC": 6.162871642559987,
        "EC": 1.988838862319425e-11,
        "NO": 0.13154598641111712,
        "SB": 1.5768783011372043,
        "SM": 2.3327090578034047,
    },
    "numubar": {
        "MD": 31.039404930560103,
        "HC": 3.6271774232568497,
        "EC": 1.2504908404874962e-10,
        "NO": 0.07932357166242918,
        "SB": 0.9276169228908447,
        "SM": 1.3722895473412542,
    },
}
wmokhov = {
    "Name": "Mohkov et al. (Fermilab)",
    "tc": 182892577621.767,
    "nue": {
        "MD": 18.195884822599115,
        "HC": 8.982495802241184,
        "EC": 1.6316400235210167,
        "NO": 1.6488589759173493,
        "SB": 0.40628892516998,
        "SM": 0.6001777168878794,
    },
    "numubar": {
        "MD": 10.749190518034164,
        "HC": 5.303804243492552,
        "EC": 0.9643313285445148,
        "NO": 0.964947855314657,
        "SB": 0.23980913837832243,
        "SM": 0.35487008916862106,
    },
    "nuebar": {
        "MD": 10.726530916266308,
        "HC": 5.305839760935416,
        "EC": 0.9675137797536498,
        "NO": 0.9759483241323901,
        "SB": 0.23896391998864214,
        "SM": 0.35298668761568947,
    },
    "numu": {
        "MD": 18.147809466147024,
        "HC": 8.971637043832846,
        "EC": 1.6375645457244183,
        "NO": 1.630277823615734,
        "SB": 0.4043998924531374,
        "SM": 0.5982284002653862,
    },
}
wmucols2 = {
    "Name": "IMCC-II",
    "tc": 165593902855.0419,
    "nue": {
        "MD": 14.828737220689618,
        "HC": 7.268344410579092,
        "EC": 0.9021635445152285,  # NOTE: I changed duplicate HC to EC instead
        "NO": 7.234258229400494,
        "SB": 0.3388182233432486,
        "SM": 0.504226158453185,
    },
    "nuebar": {
        "MD": 9.039835317040788,
        "HC": 4.411161954111631,
        "EC": 0.551511030202595,
        "NO": 4.4208021588850475,
        "SB": 0.2072642492917675,
        "SM": 0.30793775563130443,
    },
    "numu": {
        "MD": 14.836400495678491,
        "HC": 7.241186746593437,
        "EC": 0.9108501350191356,
        "NO": 7.265860179261139,
        "SB": 0.34134355096691377,
        "SM": 0.507198602315641,
    },
    "numubar": {
        "MD": 9.002308594040738,
        "HC": 4.4187549930826115,
        "EC": 0.5465999503934863,
        "NO": 4.400229459832773,
        "SB": 0.20636824645559132,
        "SM": 0.30783879421605,
    },
}

weights = {"muTs": wmuTs, "mokhov": wmokhov, "mucols2": wmucols2}

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


#density in g/cm**3; atomic mass in g/mol
Si = subs(2.329, 28.0855, 28, 14)
WSi2 = subs(9.3, 240.01, 240, 102)
Fe = subs(7.874, 55.845,56, 26)
Al = subs(2.7,26.981539, 27, 13)
W = subs(19.3, 183.84, 184, 74)
Cu = subs(8.96, 63.546, 64, 29)
PS = subs(1.05, 104.1, 104, 56)


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
    "beam_dtheta": 1.70e-4,
    "circular": False,
    "Nmu": (1.4e10 * 3.6/2.3 * (1-1/np.e) * 40) * (50) * (365.25 * 24 * 3600 * 0.7), #3.6e-9 / (1.6e-19) * (1 - 1 / np.e) * 40 * 365.25 * 24 * 3600 * 50 * 0.7,
    "syr": 365.25 * 24 * 3600 * 0.7,
    "bunch": 40,
    "finj": 50,
    "Lss": 75,
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
    "Lss": 100,
}
scd_cern = {
    "name": "SCD CERN",
    "beam_p0": 5e3,
    "pmax": 15e3,
    "pmin": 0,
    "C": 866700,
    "beam_dpop": 1e-3,
    "beam_dtheta": COMMON_dTHETA,  # 5.88e-4,
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
    "beam_dtheta": COMMON_dTHETA,  # 5.93e-4,
    "circular": False,
    "Nmu": (2e12 * 2) * (15) * (200 * 24 *3600),# 1.28e10 * 2730 * 139 / 365 * 3.154e7,
    "syr": 1.2e7,
    "bunch": 1,
    "finj": 15,
    "Lss": 50,
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
    "tc": 81264561993.55392,
    "nue":{'MD': 52.75018261623187, 'HC': 6.162871642559987, 'EC': 1.988838862319425e-11, 'NO': 0.13154598641111712, 'SB': 1.5768783011372043, 'SM': 2.3327090578034047},
    "numubar": {'MD': 31.039404930560103, 'HC': 3.6271774232568497, 'EC': 1.2504908404874962e-10, 'NO': 0.07932357166242918, 'SB': 0.9276169228908447, 'SM': 1.3722895473412542},
}
wmokhov = {
    "Name": "Mohkov et al. (Fermilab)",
    "tc": 182892577621.767,
    "nue": {'MD': 18.195884822599115, 'HC': 8.982495802241184, 'EC': 1.6316400235210167, 'NO': 1.6488589759173493, 'SB': 0.40628892516998, 'SM': 0.6001777168878794},
    "numubar": {'MD': 10.749190518034164, 'HC': 5.303804243492552, 'EC': 0.9643313285445148, 'NO': 0.964947855314657, 'SB': 0.23980913837832243, 'SM': 0.35487008916862106},
    "nuebar": {'MD': 10.726530916266308, 'HC': 5.305839760935416, 'EC': 0.9675137797536498, 'NO': 0.9759483241323901, 'SB': 0.23896391998864214, 'SM': 0.35298668761568947},
    "numu": {'MD': 18.147809466147024, 'HC': 8.971637043832846, 'EC': 1.6375645457244183, 'NO': 1.630277823615734, 'SB': 0.4043998924531374, 'SM': 0.5982284002653862},
}
wmucols2 = {
    "Name": "IMCC-II",
    "tc": 165593902855.0419,
    "nue": {'MD': 14.828737220689618, 'HC': 7.268344410579092, 'HC': 0.9021635445152285, 'NO': 7.234258229400494, 'SB': 0.3388182233432486, 'SM': 0.504226158453185},
    "nuebar":{'MD': 9.039835317040788, 'HC': 4.411161954111631, 'EC': 0.551511030202595, 'NO': 4.4208021588850475, 'SB': 0.2072642492917675, 'SM': 0.30793775563130443},
    "numu":{'MD': 14.836400495678491, 'HC': 7.241186746593437, 'EC': 0.9108501350191356, 'NO': 7.265860179261139, 'SB': 0.34134355096691377, 'SM': 0.507198602315641},
    "numubar":{'MD': 9.002308594040738, 'HC': 4.4187549930826115, 'EC': 0.5465999503934863, 'NO': 4.400229459832773, 'SB': 0.20636824645559132, 'SM': 0.30783879421605},
}

weights = {"muTs": wmuTs, "mokhov": wmokhov, "mucols2": wmucols2}

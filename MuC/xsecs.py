import numpy as np
from scipy.interpolate import interp1d

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

from DarkNews import const


"""
 Cross sections taken from Alfonso's paper

 * we will eventually assume isoscalar targets, so averaging p and n
 * everything is in cm^2
 * coherent channels are calculated for C or O, and then divided by number of total nucleons 
    (we approximate the xsec on other targets by ratios of Z_C^2/Z_X^2, with C carbon and X another nucleus).
"""

kwargs_interp = {"kind": "linear", "fill_value": 0, "bounds_error": False}

data_p = np.genfromtxt(
    files("MuC.xsec_data").joinpath("Alfonso_data_proton.dat").open()
).T
data_n = np.genfromtxt(
    files("MuC.xsec_data").joinpath("Alfonso_data_neutron.dat").open()
).T
Enu_data = data_p[0]

sigma_lightquark_CC_nu = interp1d(
    Enu_data, (data_n[1] + data_p[1]) / 2 * 1e-36, **kwargs_interp
)
sigma_lightquark_CC_nubar = interp1d(
    data_p[0], (data_n[6] + data_p[6]) / 2 * 1e-36, **kwargs_interp
)
sigma_charm_CC_nu = interp1d(
    Enu_data,
    (data_n[1] * data_n[3] + data_p[1] * data_p[3]) / 2 * 1e-36,
    **kwargs_interp,
)
sigma_charm_CC_nubar = interp1d(
    Enu_data,
    (data_n[6] * data_n[8] + data_p[6] * data_p[8]) / 2 * 1e-36,
    **kwargs_interp,
)

# Neutral current -- NOTE: approximate
sigma_NC_nu = lambda x: 0.289 / 0.66 * sigma_lightquark_CC_nu(x)
sigma_NC_nubar = lambda x: 0.289 / 0.66 * sigma_lightquark_CC_nubar(x)

# Coherent pion production (calculation on Carbon)
e, ratio = np.genfromtxt(
    files("MuC.xsec_data").joinpath("tot_xsec_cohpi0_C.dat").open(), unpack=True
)
sigma_cohpi0 = interp1d(
    e,
    ratio * 1e-40 / 12,
    kind="linear",
    fill_value=ratio[-1] * 1e-40 / 12,
    bounds_error=False,
)


# Total cross sections -- neglecting nu-e and W production
# Assume (sigma_e = sigma_mu)
total_sigmanue = (
    lambda x: sigma_lightquark_CC_nu(x) + sigma_charm_CC_nu(x) + sigma_NC_nu(x)
)
total_sigmanuebar = (
    lambda x: sigma_lightquark_CC_nubar(x) + sigma_charm_CC_nubar(x) + sigma_NC_nubar(x)
)
total_sigmanumu = (
    lambda x: sigma_lightquark_CC_nu(x) + sigma_charm_CC_nu(x) + sigma_NC_nu(x)
)
total_sigmanumubar = (
    lambda x: sigma_lightquark_CC_nubar(x) + sigma_charm_CC_nubar(x) + sigma_NC_nubar(x)
)


# """Returns necessary cross neutrino cross sections from xsecs"""
# log10E, sigmaeo, sigmamuo, _, sigmaebaro, sigmamubaro, _ = np.genfromtxt(
#     files("MuC.xsec_data").joinpath("XCC.dat").open(), unpack=True
# )
# exs = 10 ** (log10E)
# exs = np.append(exs, 1e4)

# # adding values linearily up to 10 TeV
# sigmae = np.append(sigmaeo, 7.046434082801659e-1)
# sigmamu = np.append(sigmamuo, 7.046434082801659e-1)
# sigmaebar = np.append(sigmaebaro, 3.758631195493489e-1)
# sigmamubar = np.append(sigmamubaro, 3.758631195493489e-1)

# # generating interpolators
# total_sigmanue = interp1d(exs, sigmae * exs * 1e-38, bounds_error=False, fill_value=0.0)
# total_sigmanuebar = interp1d(
#     exs, sigmaebar * exs * 1e-38, bounds_error=False, fill_value=0.0
# )
# total_sigmanumu = interp1d(
#     exs, sigmamu * exs * 1e-38, bounds_error=False, fill_value=0.0
# )
# total_sigmanumubar = interp1d(
#     exs, sigmamubar * exs * 1e-38, bounds_error=False, fill_value=0.0
# )


""" 
    Neutrino-electron elastic scattering functions 
     
     * These are all calculable
     * NOTE: neglecting radiative corrections of O(alpha_EM) for now 
"""
# Neutrino-electron scattering tools
C_LL_dic = {}
C_LL_dic["e"] = 0.7276
C_LL_dic["m"] = -0.2730
C_LL_dic["t"] = -0.2730


def nue_elastic_sigma(Enu, nuf, anu=False, C_LR=0.2334, m=const.m_e):
    """
    Calculate the total cross section ignoring terms proportional to alpha_EM.

    Parameters:
    E_nu : float
        Neutrino energy.
    C_LL : float
        Coupling constant for LL component.
    C_LR : float
        Coupling constant for LR component.
    m : float
        Mass of the lepton in GeV.

    Returns:
    tuple: Total cross sections for (neutrino-electron -> neutrino-lepton-electron) and (antineutrino-electron -> antineutrino-lepton-electron).
    """

    C_LL = C_LL_dic[nuf]
    s = 2 * m * Enu  # + m**2
    prefactor = const.Gf**2 * s / np.pi

    y_max = 2 * Enu / (2 * Enu + m)

    if not anu:
        term1 = (C_LL**2) * y_max
        term2 = (C_LR**2) * (y_max - y_max**2 + y_max**3 / 3)
        term3 = (-C_LL * C_LR * m / (2 * Enu)) * y_max**2

        sigma = prefactor * (term1 + term2 + term3)
    else:

        term1 = (C_LR**2) * y_max
        term2 = (C_LL**2) * (y_max - y_max**2 + y_max**3 / 3)
        term3 = (-C_LL * C_LR * m / (2 * Enu)) * y_max**2
        sigma = prefactor * (term1 + term2 + term3)

    return sigma * const.invGeV2_to_cm2


def nue_elastic_dsigma_dy(y, Enu, nuf, C_LR=0.2334, anu=True, m=const.m_e):
    """
    Calculate the differential cross section ignoring terms proportional to alpha_EM.

    Parameters:
    C_LL : float
        Coupling constant for LL component.
    C_LR : float
        Coupling constant for LR component.
    E_nu : float
        Neutrino energy.
    y : float
        Fractional energy transfer.
    m : float
        Mass of the lepton.

    Returns:
    tuple: Differential cross sections for (neutrino-electron -> neutrino-lepton-electron) and (antineutrino-electron -> antineutrino-lepton-electron).
    """
    C_LL = C_LL_dic[nuf]
    s = 2 * m * Enu + m**2

    if not anu:
        # Differential cross section for neutrino scattering
        dsigma_dy = (const.Gf**2 * s / np.pi) * (
            (C_LL**2) + (C_LR**2) * (1 - y) ** 2 - (C_LL * C_LR * m * y) / Enu
        )
    else:
        # Differential cross section for antineutrino scattering
        dsigma_dy = (const.Gf**2 * s / np.pi) * (
            (C_LR**2) + (C_LL**2) * (1 - y) ** 2 - (C_LL * C_LR * m * y) / Enu
        )

    return dsigma_dy


""" 
    Neutrino-nucleus trident scattering functions 
     
     * Calculated by B. Zhou and J. Beacom
     * NOTE: neglecting radiative corrections of O(alpha_EM) for now 
"""

flavors = ["e", "m", "l"]
targets = ["O16"]
colors = iter(
    ["lightblue", "dodgerblue", "blue", "lightgreen", "darkgreen", "darkgreen"]
)


def get_trident_xsec(nui, nuf, minus_lep, plus_lep, target):
    """
    Get the trident cross section for a given channel.

    Parameters:
    nui : str
        Initial neutrino flavor.
    nuf : str
        Final neutrino flavor.
    minus_lep : str
        Negative charged lepton flavor.
    plus_lep : str
        Positive lepton flavor.
    target : str
        Target nucleus.

    NOTE: neutrino and antineutrino channels (CP conjugates) have the same xsec
    E.g.,
        CP (nui --> nuf + minus_lep + plus_lep) = nuibar --> nufbar + plus_lep + minus_lep

    Returns:
    tuple: Energy and cross section for the given channel.
    """
    try:
        xsec = np.loadtxt(
            f"rare_xsecs/trident_production/v{nui}{target}TOv{nuf}{minus_lep}{plus_lep.capitalize()}X_tot.dat"
        ).T
    except FileNotFoundError:
        try:
            xsec = np.loadtxt(
                f"rare_xsecs/trident_production/v{nui}{target}TO{minus_lep}v{nuf}{plus_lep.capitalize()}X_tot.dat"
            ).T
        except FileNotFoundError:
            return None

    return interp1d(xsec[0], xsec[1], **kwargs_interp)

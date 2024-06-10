# useful densities and cross sections and stuff

import numpy as np
from scipy import interpolate

def cs_interp():
    log10E,sigmaeo,sigmamuo,_,sigmaebaro,sigmamubaro,_ = np.genfromtxt('xsecs/XCC.dat',unpack=True)
    exs = 10**(log10E)
    exs=np.append(exs, 1e4)
    sigmae = np.append(sigmaeo, 7.046434082801659e-1)
    sigmamu = np.append(sigmamuo, 7.046434082801659e-1)
    sigmaebar = np.append(sigmaebaro, 3.758631195493489e-1)
    sigmamubar = np.append(sigmamubaro, 3.758631195493489e-1)

    sigmanue = interpolate.interp1d(exs,sigmae*exs*1e-38,bounds_error=False,fill_value=0.0)
    sigmanuebar = interpolate.interp1d(exs,sigmaebar*exs*1e-38,bounds_error=False,fill_value=0.0)
    sigmanumu = interpolate.interp1d(exs,sigmamu*exs*1e-38,bounds_error=False,fill_value=0.0)
    sigmanumubar = interpolate.interp1d(exs,sigmamubar*exs*1e-38,bounds_error=False,fill_value=0.0)
    
    return sigmanue, sigmanuebar, sigmanumu, sigmanumubar

#density in g/cm**3
#atomic mass in g/mol
materials = {"densities": {"Si":2.329, "WSi2": 9.3, "Fe": 7.874, "Al":2.7, "W": 19.28}, "am": {"Si": 28.0855, "WSi2": 240.01, "Fe": 55.845, "Al":26.981539, "W": 183.84}, "Z":{"Si": 28, "WSi2": 240, "Fe": 56, "Al":27, "W": 184}}
AVOGADRO = 6.02214e23

#parameters
mutristan_small = {"beam_p0": 1e3, "pmax":3e3, "pmin":0,"Racc": 3e5/2/np.pi, "beam_dpop": 0.1, "circular": True}
mutristan_large = {"beam_p0": 3e3, "pmax":9e3, "pmin":0,"Racc": 9e5/2/np.pi, "beam_dpop": 0.1, "circular": True}

parameters = {"mutristan_small": mutristan_small, "mutristan_large": mutristan_large}

def get_nd(element):
    return AVOGADRO * materials["densities"][element]/ materials["am"][element] * materials["Z"][element]
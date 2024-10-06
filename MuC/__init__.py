import importlib
from types import ModuleType

from MuC import integrands
from MuC import mudecay_tools
from MuC import detector_tools
from MuC import collider_tools
from MuC import detgeo

# Relevant dictionaries for treating sim demands.
anti_neutrinos = ["nuebar", "numubar"]
neutrinos = ["numu", "nue"]
part_names = {"nue": "ν_e", "nuebar": "anti ν_e", "numu": "ν_μ", "numubar": "anti ν_μ"}
partn_names = {"12": "ν_e", "-12": "anti ν_e", "14": "ν_μ", "-14": "anti ν_μ"}


directions = ["left", "left", "right", "right"]
compsto2 = {
    "muon_detector": "MD",
    "solenoid_borders": "SB",
    "solenoid_mid": "SM",
    "hcal": "HC",
    "ecal": "EC",
    "nozzles": "NO",
}
pdg2names = {
    "12": "nue",
    "-12": "nuebar",
    "14": "numu",
    "-14": "numubar",
    "16": "nutau",
    "-16": "nutaubar",
}


def safe_value(key, value):
    if callable(value) or "__" in key or isinstance(value, ModuleType):
        return None
    else:
        return value


# This turns our detector modules into classes to make life easier
class Geom:
    def __init__(self, module):
        for name in dir(module):
            key, value = name, getattr(module, name)
            setattr(self, key, safe_value(key, value))


def turn_module_to_class(module):
    return Geom(module)


det_v1 = turn_module_to_class(importlib.import_module("MuC.detector_geometries.det_v1"))
det_v2 = turn_module_to_class(importlib.import_module("MuC.detector_geometries.det_v2"))
block_test = turn_module_to_class(
    importlib.import_module("MuC.detector_geometries.block_test")
)
nunulum = turn_module_to_class(
    importlib.import_module("MuC.detector_geometries.nunulum")
)

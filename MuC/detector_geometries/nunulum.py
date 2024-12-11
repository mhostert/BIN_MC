import MuC.detector_tools as det

# initializing detector components
xyplane = det.Component(0)
outside = det.Component(0)


# non-detector
zero = det.BarrelCap(outside, 0, [], 0.0, 0.0, 0.0)  # needs only the special = "end"

# initializers
one = det.InitialFaces(xyplane, 1, [], 0, 1.5)

# All objects
OBJECTS = [zero, one]

# Detector specifications, required iterations, and face_dict
zbeginning = 0
rmax = 1.5
rbp = 2.2
zending = -1 * zbeginning
INITIALIZERS = [one]  # should not contain minus_one or zero or minus_two
DECAYER = []
OUTSIDE = []
iterations = 2
facedict = {}

name = "Nu Nu Luminosity"

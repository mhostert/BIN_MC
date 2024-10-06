from MuC.detector_tools import face, decayer, cap, initializer, barrel, conic
import MuC.collider_tools as ud
import copy

# initializing detector components
xyplane = face(0)
outside = face(0)


# non-detector
zero = cap(outside, 0, [], 0.0, 0.0, 0.0)  # needs only the special = "end"

# initializers
one = initializer(xyplane, 1, [], 0, 1.5)

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

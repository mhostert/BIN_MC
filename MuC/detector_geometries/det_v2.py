import MuC.detector_tools as det
import copy

# initializing detector components
muon_detector = det.Component(det.Fe.N)
solenoid_borders = det.Component(det.Fe.N)
solenoid_mid = det.Component(det.Al.N)
hcal = det.Component(det.hcal_CLICdet.N)  # Fe absorber, middle parts comple
ecal = det.Component(
    det.ecal_CLICdet.N
)  # W absrober plates, silicon detector, have to weigh them appropriately
vacuum_in_between = det.Component(det.Air.N)  # this includes the very thin tracker
nozzle = det.Component(det.W.N)
beampipe = det.Component(0)
outside = det.Component(0)

###### decayers
# decays outside beampipe, if any
minus_two = det.MuonContainer(parent=nozzle, id=-2, next_ids=[13, 50], last=50)
# decays in beampipe
minus_one = det.MuonContainer(
    parent=beampipe, id=-1, next_ids=[15, 16, 17, 18], last=18
)

###### non-detector
zero = det.BarrelCap(outside, 0, [], 0.0, 0.0, 0.0)  # needs only the special = "end"

###### initializers
one = det.InitialFaces(
    parent=muon_detector, id=1, next_ids=[4, 6, 88], rsmall=78.2, rbig=645
)
two = det.InitialFaces(
    parent=nozzle, id=2, next_ids=[8, 9, 10, 11, 12, 13], rsmall=2.2, rbig=78.2
)
three = det.InitialFaces(
    parent=beampipe, id=3, next_ids=[15, 16, 17, 18], rsmall=0.0, rbig=2.2
)

###### detector components

four = det.Barrel(outside, 4, [], 645.0, -563.8, 563.8)
five = det.BarrelCap(outside, 5, [], 563.8, 78.2, 645.0)
six = det.BarrelCap(
    vacuum_in_between, 6, [19, 20, 21, 22, 26, 89, 90, 91, 92], -417.9, 57.5, 446.1
)
seven = det.Barrel(vacuum_in_between, 7, [22, 26, 31, 34], 446.1, -417.9, 417.9)

eightyeight = det.BarrelCap(muon_detector, 88, [4, 7, 93], -417.9, 446.1, 645.0)
eightynine = det.BarrelCap(
    vacuum_in_between, 89, [23, 28, 30, 31], -412.9, 324.6, 348.3
)
ninety = det.BarrelCap(vacuum_in_between, 90, [24, 31, 32], -412.9, 352.3, 364.9)
ninetyone = det.BarrelCap(vacuum_in_between, 91, [25, 31, 33], -412.9, 399.3, 425.0)
ninetytwo = det.BarrelCap(vacuum_in_between, 92, [26, 31, 34], -412.9, 429.0, 446.1)
ninetythree = det.BarrelCap(muon_detector, 93, [4, 5], 417.9, 446.1, 645.0)

eight = det.Conic(
    muon_detector,
    8,
    [4, 6],
    (78.2 - 57.5) / (563.8 - 417.9),
    -563.8,
    -417.9,
    57.5,
    "towards",
)
fourteen = copy.deepcopy(eight)
fourteen.id = 14
fourteen.density = nozzle.density
fourteen.next_ids = [8, 9, 10, 11, 12, 13]
nine = det.Conic(
    vacuum_in_between,
    9,
    [19],
    (57.5 - 56.8) / (417.9 - 412.9),
    -417.9,
    -412.9,
    56.8,
    "towards",
)
ten = det.Conic(
    hcal,
    10,
    [35, 36, 37],
    (56.8 - 33.9) / (412.9 - 250.9),
    -412.9,
    -250.9,
    33.9,
    "towards",
)
eleven = det.Conic(
    ecal, 11, [38], (33.9 - 31) / (250.9 - 230.7), -250.9, -230.7, 31.0, "towards"
)
twelve = det.Conic(
    vacuum_in_between,
    12,
    [41, 42, 43, 44],
    (31 - 2.2) / (230.7 - 6.5),
    -230.7,
    -6.5,
    2.2,
    "towards",
)
thirteen = det.Barrel(beampipe, 13, [15, 16, 17, 18], 2.2, -563.8, -6.5)

fifteen = det.Barrel(nozzle, 15, [8, 9, 10, 11, 12], 2.2, -563.8, -6.5)
sixteen = det.Barrel(vacuum_in_between, 16, [41, 42, 43], 2.2, -6.5, 6.5)
seventeen = det.Barrel(nozzle, 17, [45, 46, 47, 48, 49, 50], 2.2, 6.5, 563.8)
eighteen = det.BarrelCap(outside, 18, [], 563.8, 0, 2.2)

nineteen = det.BarrelCap(hcal, 19, [35, 36, 37, 52, 53], -412.9, 56.8, 324.6)
twenty = det.BarrelCap(solenoid_borders, 20, [54, 55, 56], -412.9, 348.3, 352.3)
twentyone = det.BarrelCap(solenoid_mid, 21, [57, 58, 59], -412.9, 364.9, 399.3)
twentytwo = det.BarrelCap(solenoid_borders, 22, [60, 61, 62], -412.9, 425.0, 429.0)
twentythree = det.Barrel(solenoid_borders, 23, [54, 56], 348.3, -412.9, 412.9)
twentyfour = det.Barrel(solenoid_mid, 24, [57, 59], 364.9, -412.9, 412.9)
twentyfive = det.Barrel(solenoid_borders, 25, [60, 62], 425.0, -412.9, 412.9)
twentysix = det.Barrel(muon_detector, 26, [4, 93], 446.1, -417.9, 417.9)

twentyseven = copy.deepcopy(nine)
twentyseven.id = 27
twentyseven.density = nozzle.density
twentyseven.next_ids = [9, 10, 11, 12, 13]
twentyeight = det.Barrel(hcal, 28, [35, 36, 37, 52, 53], 324.6, -412.9, -235.4)
twentynine = det.BarrelCap(hcal, 29, [63, 64, 65], -221, 174, 333)
thirty = det.Barrel(hcal, 30, [63, 64, 65], 333, -221, 221)
thirtyone = det.BarrelCap(muon_detector, 31, [5, 66], 417.9, 57.5, 446.1)
thirtytwo = det.Barrel(solenoid_borders, 32, [54, 55, 56], 352.3, -412.9, 412.9)
thirtythree = det.Barrel(solenoid_mid, 33, [57, 58, 59], 399.3, -412.9, 412.9)
thirtyfour = det.Barrel(solenoid_borders, 34, [60, 61, 62], 429.0, -412.9, 412.9)

thirtyfive = det.BarrelCap(ecal, 35, [38, 39, 40, 67], -250.9, 33.9, 170.0)
thirtysix = det.BarrelCap(vacuum_in_between, 36, [29, 71], -235.4, 170.0, 324.6)
thirtyseven = det.Barrel(vacuum_in_between, 37, [23, 29], 324.6, -412.9, -235.4)
thirtyeight = det.BarrelCap(
    vacuum_in_between, 38, [41, 42, 43, 44, 69, 70, 71, 72, 87], -230.7, 31.0, 170.0
)
thirtynine = det.Barrel(hcal, 39, [36], 170.0, -250.9, -235.4)
forty = det.Barrel(vacuum_in_between, 40, [23, 29, 69, 71], 170.0, -235.4, -230.7)

fortyone = det.Barrel(ecal, 41, [73, 74], 150, -221, 221)
fortytwo = det.BarrelCap(ecal, 42, [75, 76, 77, 78], 230.7, 31, 170)
fortythree = det.Conic(
    nozzle,
    43,
    [45, 46, 47, 48, 49, 50, 51],
    (31 - 2.2) / (230.7 - 6.5),
    6.5,
    230.7,
    2.2,
    "away",
)
fortyfour = det.Barrel(beampipe, 44, [16, 17, 18], 2.2, -6.5, 6.5)
fortyfive = copy.deepcopy(fortythree)
fortyfive.id = 45
fortyfive.density = vacuum_in_between.density
fortyfive.next_ids = [41, 42]

seventyeight = det.Conic(
    nozzle,
    78,
    [46, 47, 49, 50, 51],
    (33.9 - 31) / (250.9 - 230.7),
    230.7,
    250.9,
    31.0,
    "away",
)
fortysix = copy.deepcopy(seventyeight)
fortysix.id = 46
fortysix.density = ecal.density
fortysix.next_ids = [77]
fortyseven = det.Conic(
    hcal, 47, [79, 80], (56.8 - 33.9) / (412.9 - 250.9), 250.9, 412.9, 33.9, "away"
)
fortyeight = det.Conic(
    vacuum_in_between,
    48,
    [31],
    (57.5 - 56.8) / (417.9 - 412.9),
    412.9,
    417.9,
    56.8,
    "away",
)
fortynine = det.Conic(
    muon_detector, 49, [5], (78.2 - 57.5) / (563.8 - 417.9), 417.9, 563.8, 57.5, "away"
)
fifty = det.BarrelCap(outside, 50, [], 563.8, 2.2, 78.2)
fiftyone = det.Barrel(beampipe, 51, [17, 18], 2.2, 6.5, 563.8)

fiftytwo = det.Barrel(ecal, 52, [38, 67], 170.0, -250.9, -235.4)
fiftythree = copy.deepcopy(ten)
fiftythree.id = 53
fiftythree.density = nozzle.density
fiftythree.next_ids = [10, 11, 12, 13]
fiftyfour = det.Barrel(vacuum_in_between, 54, [24, 31], 352.3, -412.9, 412.9)
fiftyfive = det.Barrel(
    vacuum_in_between, 55, [23, 28, 29, 30, 31, 81], 348.3, -412.9, 412.9
)
fiftysix = det.BarrelCap(vacuum_in_between, 56, [31], 412.9, 348.3, 352.3)
fiftyseven = det.Barrel(vacuum_in_between, 57, [25, 31], 399.3, -412.9, 412.9)

fiftyeight = det.Barrel(vacuum_in_between, 58, [24, 31, 32], 364.9, -412.9, 412.9)
fiftynine = det.BarrelCap(vacuum_in_between, 59, [31], 412.9, 364.9, 399.3)
sixty = det.Barrel(vacuum_in_between, 60, [26, 31], 429.0, -412.9, 412.9)
sixtyone = det.Barrel(vacuum_in_between, 61, [25, 31, 33], 425.0, -412.9, 412.9)
sixtytwo = det.BarrelCap(vacuum_in_between, 62, [31], 412.9, 425.0, 429.0)

sixtythree = det.Barrel(vacuum_in_between, 63, [23, 31], 333.0, -221.0, 221.0)
sixtyfour = det.BarrelCap(vacuum_in_between, 64, [23, 31, 72], 221.0, 174.0, 333.0)
sixtyfive = det.Barrel(vacuum_in_between, 65, [70, 72, 83], 174.0, -221.0, 221.0)
sixtysix = copy.deepcopy(fortynine)
sixtysix.id = 66
sixtysix.density = nozzle.density
sixtysix.next_ids = [49, 50, 51]

sixtyseven = copy.deepcopy(eleven)
sixtyseven.id = 67
sixtyseven.density = nozzle.density
sixtyseven.next_ids = [11, 12, 13]
sixtyeight = copy.deepcopy(forty)
sixtyeight.id = 68
sixtyeight.density = ecal.density
sixtyeight.next_ids = [35]
sixtynine = det.BarrelCap(ecal, 69, [73, 74, 84], -221.0, 150.0, 170.2)
seventy = det.Barrel(ecal, 70, [73, 74, 84], 170.2, -221.0, 221.0)

seventyone = copy.deepcopy(sixtyfive)
seventyone.id = 71
seventyone.density = hcal.density
seventyone.next_ids = [63, 64]
seventytwo = det.BarrelCap(hcal, 72, [79, 80, 85], 235.4, 170.0, 324.6)
seventythree = copy.deepcopy(seventy)
seventythree.id = 73
seventythree.density = vacuum_in_between.density
seventythree.next_ids = [71, 72]

seventyfour = det.BarrelCap(vacuum_in_between, 74, [42, 71, 72], 221.0, 150.0, 170.2)
seventyfive = det.Barrel(vacuum_in_between, 75, [72], 170.0, 230.7, 235.4)
seventysix = det.Barrel(hcal, 76, [79, 80], 170.0, 235.4, 250.9)
seventyseven = det.BarrelCap(hcal, 77, [79, 80, 85], 250.9, 33.9, 170.0)
seventynine = det.BarrelCap(vacuum_in_between, 79, [31, 82], 412.9, 56.8, 324.6)
eighty = det.Barrel(vacuum_in_between, 80, [23, 31], 324.6, 235.4, 412.9)

eightyone = copy.deepcopy(eighty)
eightyone.id = 81
eightyone.density = hcal.density
eightyone.next_ids = [79, 80, 85]
eightytwo = copy.deepcopy(fortyeight)
eightytwo.id = 82
eightytwo.density = nozzle.density
eightytwo.next_ids = [48, 49, 50, 51]
eightythree = copy.deepcopy(seventyfive)
eightythree.id = 83
eightythree.density = ecal.density
eightythree.next_ids = [77]

eightyfour = copy.deepcopy(fortyone)
eightyfour.id = 84
eightyfour.density = vacuum_in_between.density
eightyfour.next_ids = [41, 42, 43, 44, 87]
eightyfive = copy.deepcopy(fortyseven)
eightyfive.id = 85
eightyfive.density = nozzle.density
eightyfive.next_ids = [47, 48, 49, 50, 51]

eightysix = copy.deepcopy(seventysix)
eightysix.id = 86
eightysix.density = ecal.density
eightysix.next_ids = [77]
eightyseven = copy.deepcopy(twelve)
eightyseven.id = 87
eightyseven.density = nozzle.density
eightyseven.next_ids = [12, 13]

# All objects
OBJECTS = [
    zero,
    one,
    two,
    three,
    four,
    five,
    six,
    seven,
    eight,
    nine,
    ten,
    eleven,
    twelve,
    thirteen,
    fourteen,
    fifteen,
    sixteen,
    seventeen,
    eighteen,
    nineteen,
    twenty,
    twentyone,
    twentytwo,
    twentythree,
    twentyfour,
    twentyfive,
    twentysix,
    twentyseven,
    twentyeight,
    twentynine,
    thirty,
    thirtyone,
    thirtytwo,
    thirtythree,
    thirtyfour,
    thirtyfive,
    thirtysix,
    thirtyseven,
    thirtyeight,
    thirtynine,
    forty,
    fortyone,
    fortytwo,
    fortythree,
    fortyfour,
    fortyfive,
    fortysix,
    fortyseven,
    fortyeight,
    fortynine,
    fifty,
    fiftyone,
    fiftytwo,
    fiftythree,
    fiftyfour,
    fiftyfive,
    fiftysix,
    fiftyseven,
    fiftyeight,
    fiftynine,
    sixty,
    sixtyone,
    sixtytwo,
    sixtythree,
    sixtyfour,
    sixtyfive,
    sixtysix,
    sixtyseven,
    sixtyeight,
    sixtynine,
    seventy,
    seventyone,
    seventytwo,
    seventythree,
    seventyfour,
    seventyfive,
    seventysix,
    seventyseven,
    seventyeight,
    seventynine,
    eighty,
    eightyone,
    eightytwo,
    eightythree,
    eightyfour,
    eightyfive,
    eightysix,
    eightyseven,
    eightyeight,
    eightynine,
    ninety,
    ninetyone,
    ninetytwo,
    ninetythree,
]

# Detector specifications, required iterations, and face_dict
zbeginning = -563.8  # cm
rmax = 645.0  # cm
beam_pipe_radius = 2.2  # cm
zending = -1 * zbeginning  # cm
INITIALIZERS = [three, two, one]  # should not contain minus_one or zero or minus_two
DECAYER = [minus_one, minus_two]
OUTSIDE = [minus_one, zero, four, five, fifty, eighteen]
iterations = 35
facedict = {
    "muon_detector": [1, 31, 8, 49, 88, 26, 93],
    "solenoid_borders": [20, 22, 23, 25, 32, 34],
    "solenoid_mid": [21, 24, 33],
    "hcal": [10, 19, 28, 29, 30, 39, 47, 71, 72, 76, 77, 81],
    "ecal": [11, 35, 41, 42, 46, 52, 68, 69, 70, 83, 86],
    "nozzles": [-2, 2, 14, 15, 17, 27, 43, 78, 53, 66, 67, 82, 85, 87],
}

"""facedict = {'muon_detector_ec': [1, 31, 8, 49, 93],
            'muon_detector_barrel': [26, 88],
            'solenoid_borders': [20,22,23,25,32,34],
            'solenoid_mid':[21, 24, 33],
            'hcal_ec': [10,19,28,39,47,72,76,77,81],
            'hcal_barrel': [29,30,71],
            'ecal_ec': [11,35,42,46,52,68,83,86],
            'ecal_barrel': [41,69,70],
            'nozzles':[-2,2,14,15,17,27,43,78,53,66,67,82,85,87]}"""

name = "Det v2"

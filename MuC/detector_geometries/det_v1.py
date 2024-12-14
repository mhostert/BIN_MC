# approximate muon detector (v1)
from MuC import detector_tools as det

# these densities will be number densities
muon_detector_hcal = det.Fe
solenoid = det.Al
ecal = det.WSi2  # should actually use SiW, not WSi2
tracker = det.Si
cone = det.W
beampipe = det.vacuum
outside = det.vacuum

# caps have: id, next_ids, zpos, rbeg, rend

# barrels have: id, next_ids, rpos, zbeg, zend

# conic have: id, next_ids, tan_theta, zbeg, zend, rsmall, direction

# initializers have: id, next_ids, rsmall, rbig (their zpos is already defined in )

# decayers have: id, next_ids

minus_two = det.MuonContainer(cone, -2, [18, 17, 19, 20, 27, 36], 7)
minus_one = det.MuonContainer(outside, -1, [4, 5, 6, 7], 7)

zero = det.BarrelCap(outside, 0, [], 0, 0, 0)  # needs only the special = "end"

one = det.InitialFaces(
    muon_detector_hcal, 1, [8, 9, 10, 11, 12, 13, 14, 15, 16], 58, 645
)
two = det.InitialFaces(cone, 2, [17, 18, 19, 20], 3, 58)
three = det.InitialFaces(beampipe, 3, [4, 5, 6, 7], 0, 3)

four = det.Barrel(
    cone, 4, [17, 19, 20], 3, -564, -28
)  # four five six could be made easier
five = det.Barrel(tracker, 5, [21, 22, 23], 3, -28, 28)
six = det.Barrel(cone, 6, [24, 25, 26, 27], 3, 28, 564)
seven = det.BarrelCap(outside, 7, [], 564, 0, 3)

eight = det.Conic(cone, 8, [17, 18, 19, 20], 32 / 313, -564, -251, 26, "towards")
nine = det.BarrelCap(solenoid, 9, [28, 29, 30], -418, 348, 446)
ten = det.Barrel(solenoid, 10, [28, 30], 348, -418, 418)
eleven = det.Barrel(solenoid, 11, [28, 29, 30], 446, -418, 418)
twelve = det.BarrelCap(ecal, 12, [31, 32, 33, 34, 35], -251, 26, 170)
thirteen = det.Barrel(ecal, 13, [31, 32, 33, 34, 35, 39], 170, -251, 251)
fourteen = det.Conic(cone, 14, [26, 27, 36], 32 / 313, 251, 564, 26, "away")
fifteen = det.BarrelCap(outside, 15, [], 564, 58, 645)
sixteen = det.Barrel(outside, 16, [], 645, -564, 564)

seventeen = det.Conic(
    muon_detector_hcal, 17, [9, 10, 12, 15, 16], 32 / 313, -564, -251, 26, "towards"
)  # seventeen, nineteen, and twenty could be optimized by using a check function instead of three finding intersection with cones, but only works with same slope
eighteen = det.Barrel(beampipe, 18, [4, 5, 6, 7], 3, -564, -28)
nineteen = det.Conic(ecal, 19, [32, 33], 2 / 20, -251, -231, 24, "towards")
twenty = det.Conic(tracker, 20, [21, 22, 23, 37], 21 / 203, -231, -28, 3, "towards")

twentyone = det.Barrel(ecal, 21, [33, 35], 150, -231, 231)
twentytwo = det.BarrelCap(ecal, 22, [33, 35, 39], 231, 24, 150)
twentythree = det.Conic(cone, 23, [24, 25, 26, 27, 36], 21 / 203, 28, 231, 3, "away")

twentyfour = det.Conic(tracker, 24, [21, 22], 21 / 203, 28, 231, 3, "away")
twentyfive = det.Conic(ecal, 25, [33, 35], 2 / 20, 231, 251, 24, "away")
twentysix = det.Conic(
    muon_detector_hcal, 26, [10, 15, 16], 32 / 313, 251, 564, 26, "away"
)
twentyseven = det.BarrelCap(outside, 27, [], 564, 3, 58)

twentyeight = det.Barrel(muon_detector_hcal, 28, [15, 16], 446, -418, 418)
twentynine = det.Barrel(
    muon_detector_hcal, 29, [8, 10, 12, 13, 14, 15, 16], 348, -418, 418
)
thirty = det.BarrelCap(muon_detector_hcal, 30, [14, 15, 16], 418, 348, 446)

thirtyone = det.Conic(cone, 31, [18, 19, 20], 2 / 20, -251, -231, 24, "towards")
thirtytwo = det.BarrelCap(tracker, 32, [21, 22, 23, 37, 38], -231, 24, 150)
thirtythree = det.Barrel(muon_detector_hcal, 33, [10, 15], 170, -251, 251)
thirtyfour = det.Barrel(tracker, 34, [21, 22, 23, 37, 38], 150, -231, 231)
thirtyfive = det.BarrelCap(muon_detector_hcal, 35, [10, 14, 15, 16], 251, 26, 170)

thirtysix = det.Barrel(beampipe, 36, [6, 7], 3, 28, 564)

thirtyseven = det.Barrel(beampipe, 37, [5, 6, 7], 3, -28, 28)
thirtyeight = det.Conic(cone, 38, [18, 20], 21 / 203, -231, -28, 3, "towards")

thirtynine = det.Conic(cone, 39, [25, 26, 27, 36], 2 / 20, 231, 251, 24, "away")


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
]
zbeginning = -564
rmax = 645
beam_pipe_radius = 3
zending = -1 * zbeginning
INITIALIZERS = [three, two, one]  # should not contain minus_one or zero
DECAYER = [minus_one, minus_two]
OUTSIDE = [minus_one, zero, seven, fifteen, sixteen, twentyseven]
iterations = 25

name = "Det v1"


# approximate muon detector (v1)

from helpers import face
import useful_data as ud

#these densities will be number densities
muon_detector_hcal = face(ud.get_nd("Fe"))
solenoid = face(ud.get_nd("Al"))
ecal = face(ud.get_nd("WSi2"))
tracker = face(ud.get_nd("Si"))
cone = face(ud.get_nd("W"))
beampipe = face(0)
outside = face(0)

# caps have: id, next_ids, zpos, rbeg, rend

# barrels have: id, next_ids, rpos, zbeg, zend

# conic have: id, next_ids, tan_theta, zbeg, zend, rsmall, direction

# initializers have: id, next_ids, rsmall, rbig (their zpos is already defined in )

minus_one = outside.initializer(outside,-1, [4,5,6,7], 0,0) # will not use initialzier check_in, only neighbors

zero = outside.cap(outside,0, [], 0,0,0) # needs only the special = "end"

one = muon_detector_hcal.initializer(muon_detector_hcal,1, [8,9,10,11,12,13,14,15,16], 58, 645)
two = cone.initializer(cone, 2, [17,18,19,20], 3,58)
three = beampipe.initializer(beampipe,3, [4,5,6,7], 0, 3)

four = cone.barrel(cone,4, [17,19,20], 3, -564, -28) #four five six could be made easier
five = tracker.barrel(tracker,5, [21,22,23], 3, -28, 28)
six = cone.barrel(cone,6, [24,25,26,27], 3, 28, 564)
seven = outside.cap(outside,7, [], 564, 0, 3)

eight = cone.conic(cone,8, [17,18,19,20], 32/313, -564, -251, 26, "towards")
nine = solenoid.cap(solenoid,9, [28,29,30], -418, 348, 446)
ten = solenoid.barrel(solenoid,10, [28,30], 348, -418, 418)
eleven = solenoid.barrel(solenoid,11, [28,29,30], 446, -418, 418)
twelve = ecal.cap(ecal,12, [31,32,33,34,35], -251, 26, 170)
thirteen = ecal.barrel(ecal,13, [31,32,33,34,35,39], 170, -251, 251)
fourteen = cone.conic(cone,14, [26,27,36], 32/313, 251, 564, 26, "away")
fifteen = outside.cap(outside,15, [], 564, 58,645)
sixteen = outside.barrel(outside,16, [], 645, -564, 564)

seventeen = muon_detector_hcal.conic(muon_detector_hcal,17, [9,10,12,15,16],32/313, -564, -251, 26, "towards") #seventeen, nineteen, and twenty could be optimized by using a check function instead of three finding intersection with cones, but only works with same slope
eighteen = beampipe.barrel(beampipe,18, [4,5,6,7], 3, -564, -28)
nineteen = ecal.conic(ecal,19, [32,33], 2/20 ,-251, -231, 24, "towards")
twenty = tracker.conic(tracker,20, [21,22,23, 37], 21/203, -231, -28, 3, "towards")

twentyone = ecal.barrel(ecal,21, [33,35],150, -231, 231)
twentytwo = ecal.cap(ecal,22, [33,35,39], 231, 24,150)
twentythree = cone.conic(cone,23, [24,25,26,27,36], 21/203, 28, 231, 3, "away")

twentyfour = tracker.conic(tracker,24, [21,22], 21/203, 28, 231, 3, "away")
twentyfive = ecal.conic(ecal,25, [33,35], 2/20, 231, 251, 24, "away")
twentysix = muon_detector_hcal.conic(muon_detector_hcal,26, [10,15,16], 32/313, 251, 564, 26, "away")
twentyseven = outside.cap(outside,27,[], 564, 3, 58)

twentyeight = muon_detector_hcal.barrel(muon_detector_hcal,28, [15,16], 446, -418, 418)
twentynine = muon_detector_hcal.barrel(muon_detector_hcal,29, [8,10,12,13,14,15,16], 348, -418, 418)
thirty = muon_detector_hcal.cap(muon_detector_hcal,30, [14,15,16],418, 348, 446)

thirtyone = cone.conic(cone,31, [18,19,20], 2/20, -251, -231, 24, "towards")
thirtytwo = tracker.cap(tracker,32, [21,22,23,37,38], -231, 24, 150)
thirtythree = muon_detector_hcal.barrel(muon_detector_hcal,33, [10,15], 170, -251, 251)
thirtyfour = tracker.barrel(tracker,34, [21,22,23,37,38], 150,-231,231)
thirtyfive = muon_detector_hcal.cap(muon_detector_hcal,35, [10,14,15,16], 251, 26, 170)

thirtysix = beampipe.barrel(beampipe,36, [6,7], 3, 28, 564)

thirtyseven = beampipe.barrel(beampipe,37, [5,6,7], 3, -28, 28)
thirtyeight = cone.conic(cone,38, [18,20], 21/203, -231, -28, 3, "towards")

thirtynine = cone.conic(cone,39, [25,26,27,36], 2/20,231,251,24,"away")





OBJECTS = [zero, one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen, nineteen, twenty, twentyone, twentytwo, twentythree, twentyfour, twentyfive, twentysix, twentyseven, twentyeight,twentynine,thirty, thirtyone, thirtytwo, thirtythree, thirtyfour, thirtyfive, thirtysix, thirtyseven, thirtyeight, thirtynine]
zbeginning = -564
rmax = 645
zending = -1 * zbeginning
INITIALIZERS = [one, two, three] #should not contain minus_one
DECAYER = [minus_one]
OUTSIDE = [minus_one, zero,seven, fifteen, sixteen, twentyseven]

'''objects = geom.OBJECTS
    zbeginning = geom.zbeginning
    rmax = geom.rmax
    zending = geom.zending
    initials = geom.initializers #should not contain minus_one
    decayer = geom.decayer #is the minus_one'''
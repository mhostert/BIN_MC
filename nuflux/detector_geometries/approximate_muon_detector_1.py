
# approximate muon detector (v1)

from helpers import face, decayer, initializer, conic, barrel, cap
import useful_data as ud

#these densities will be number densities
muon_detector_hcal = face(ud.Fe.N)
solenoid = face(ud.Al.N)
ecal = face(ud.WSi2.N) # should actually use SiW, not WSi2
tracker = face(ud.Si.N)
cone = face(ud.W.N)
beampipe = face(0)
outside = face(0)

# caps have: id, next_ids, zpos, rbeg, rend

# barrels have: id, next_ids, rpos, zbeg, zend

# conic have: id, next_ids, tan_theta, zbeg, zend, rsmall, direction

# initializers have: id, next_ids, rsmall, rbig (their zpos is already defined in )

# decayers have: id, next_ids

minus_two = decayer(cone, -2, [18,17,19,20, 27,36], 7)
minus_one = decayer(outside,-1, [4,5,6,7], 7)

zero = cap(outside,0, [], 0,0,0) # needs only the special = "end"

one = initializer(muon_detector_hcal,1, [8,9,10,11,12,13,14,15,16], 58, 645)
two = initializer(cone, 2, [17,18,19,20], 3,58)
three = initializer(beampipe,3, [4,5,6,7], 0, 3)

four = barrel(cone,4, [17,19,20], 3, -564, -28) #four five six could be made easier
five = barrel(tracker,5, [21,22,23], 3, -28, 28)
six = barrel(cone,6, [24,25,26,27], 3, 28, 564)
seven = cap(outside,7, [], 564, 0, 3)

eight = conic(cone,8, [17,18,19,20], 32/313, -564, -251, 26, "towards")
nine = cap(solenoid,9, [28,29,30], -418, 348, 446)
ten = barrel(solenoid,10, [28,30], 348, -418, 418)
eleven = barrel(solenoid,11, [28,29,30], 446, -418, 418)
twelve = cap(ecal,12, [31,32,33,34,35], -251, 26, 170)
thirteen = barrel(ecal,13, [31,32,33,34,35,39], 170, -251, 251)
fourteen = conic(cone,14, [26,27,36], 32/313, 251, 564, 26, "away")
fifteen = cap(outside,15, [], 564, 58,645)
sixteen = barrel(outside,16, [], 645, -564, 564)

seventeen = conic(muon_detector_hcal,17, [9,10,12,15,16],32/313, -564, -251, 26, "towards") #seventeen, nineteen, and twenty could be optimized by using a check function instead of three finding intersection with cones, but only works with same slope
eighteen = barrel(beampipe,18, [4,5,6,7], 3, -564, -28)
nineteen = conic(ecal,19, [32,33], 2/20 ,-251, -231, 24, "towards")
twenty = conic(tracker,20, [21,22,23, 37], 21/203, -231, -28, 3, "towards")

twentyone = barrel(ecal,21, [33,35],150, -231, 231)
twentytwo = cap(ecal,22, [33,35,39], 231, 24,150)
twentythree = conic(cone,23, [24,25,26,27,36], 21/203, 28, 231, 3, "away")

twentyfour = conic(tracker,24, [21,22], 21/203, 28, 231, 3, "away")
twentyfive = conic(ecal,25, [33,35], 2/20, 231, 251, 24, "away")
twentysix = conic(muon_detector_hcal,26, [10,15,16], 32/313, 251, 564, 26, "away")
twentyseven = cap(outside,27,[], 564, 3, 58)

twentyeight = barrel(muon_detector_hcal,28, [15,16], 446, -418, 418)
twentynine = barrel(muon_detector_hcal,29, [8,10,12,13,14,15,16], 348, -418, 418)
thirty = cap(muon_detector_hcal,30, [14,15,16],418, 348, 446)

thirtyone = conic(cone,31, [18,19,20], 2/20, -251, -231, 24, "towards")
thirtytwo = cap(tracker,32, [21,22,23,37,38], -231, 24, 150)
thirtythree = barrel(muon_detector_hcal,33, [10,15], 170, -251, 251)
thirtyfour = barrel(tracker,34, [21,22,23,37,38], 150,-231,231)
thirtyfive = cap(muon_detector_hcal,35, [10,14,15,16], 251, 26, 170)

thirtysix = barrel(beampipe,36, [6,7], 3, 28, 564)

thirtyseven = barrel(beampipe,37, [5,6,7], 3, -28, 28)
thirtyeight = conic(cone,38, [18,20], 21/203, -231, -28, 3, "towards")

thirtynine = conic(cone,39, [25,26,27,36], 2/20,231,251,24,"away")





OBJECTS = [zero, one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen, nineteen, twenty, twentyone, twentytwo, twentythree, twentyfour, twentyfive, twentysix, twentyseven, twentyeight,twentynine,thirty, thirtyone, thirtytwo, thirtythree, thirtyfour, thirtyfive, thirtysix, thirtyseven, thirtyeight, thirtynine]
zbeginning = -564
rmax = 645
rbp = 3
zending = -1 * zbeginning
INITIALIZERS = [three, two, one] #should not contain minus_one or zero
DECAYER = [minus_one, minus_two]
OUTSIDE = [minus_one, zero,seven, fifteen, sixteen, twentyseven]
iterations = 25
facedict = {}
TESTS = []
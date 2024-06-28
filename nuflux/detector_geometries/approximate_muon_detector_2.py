
# approximate muon detector (v2)

from helpers import face, decayer, cap, initializer, barrel, conic
import useful_data as ud
import copy

#these densities will be number densities
muon_detector = face(ud.Fe.N)
solenoid_borders = face(ud.Fe.N)
solenoid_mid = face(ud.Al.N)
hcal = face(ud.hcal_CLICdet.N) # Fe absorber, middle parts complex
ecal = face(ud.ecal_CLICdet.N) # W absrober plates, silicon detector, have to weigh them appropriately
space_between = face(0) # this includes the very thin tracker
nozzle = face(ud.W.N)
beampipe = face(0)
outside = face(0)

# caps have: id, next_ids, zpos, rbeg, rend

# barrels have: id, next_ids, rpos, zbeg, zend

# conic have: id, next_ids, tan_theta, zbeg, zend, rsmall, direction

# initializers have: id, next_ids, rsmall, rbig (their zpos is already defined in detgeo) (these are the elements of the first barrier to be in the detector)

# decayers have: id, next_ids

minus_two = decayer(nozzle, -2, [13, 50]) #decays outside beampipe, if any
minus_one = decayer(beampipe,-1, [15, 16, 17, 18]) #decays in beampipe

zero = cap(outside,0, [], 0,0,0) # needs only the special = "end"

one = initializer(muon_detector, 1, [4,5,6,7,14], 78.2, 645)
two = initializer(nozzle,2, [8,9,10,11,12,13], 2.2, 78.2)
three = initializer(beampipe, 3, [15,16,17,18], 0, 2.2)

four = barrel(outside, 4, [], 645, -563.8, 563.8)
five = cap(outside, 5, [], 563.8, 78.2, 645)
six = cap(space_between, 6, [19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34], -417.9, 57.5,  446.1)
seven = barrel(space_between, 7, [22,26,31,34], 446.1, -417.9, 417.9)



eight = conic(muon_detector, 8, [4,6], (78.2 - 57.5)/(563.8 - 417.9), -563.8, -417.9, 57.5, "towards")
fourteen = copy.deepcopy(eight)
fourteen.id = 14
fourteen.density = nozzle.density
fourteen.next_ids = [8,9,10,11,12,13]
nine = conic(space_between, 9, [19], (57.5 - 56.8)/( 417.9 - 412.9), -417.9, -412.9, 56.8, 'towards')
ten = conic(hcal, 10, [35,36,37], (56.8 - 33.9)/(412.9 - 250.9), -412.9, -250.9, 33.9, 'towards')
eleven = conic(ecal, 11, [38], (33.9 - 31)/(250.9 - 230.7), -250.9, -230.7, 31, 'towards')
twelve = conic(space_between, 12, [41,42,43,44], (31-2.2)/(230.7 - 6.5), -230.7, -6.5, 2.2, 'towards')
thirteen = barrel(beampipe, 13, [15,16,17,18], 2.2, -563.8, -6.5)

fifteen = barrel(nozzle, 15, [8,9,10,11,12], 2.2, -563.8, -6.5)
sixteen = barrel(space_between, 16, [41,42,43], 2.2, -6.5, 6.5)
seventeen = barrel(nozzle, 17, [45,46,47,48,49,50], 2.2, 6.5, 563.8)
eighteen = cap(outside, 18, [], 563.8, 0, 2.2)

nineteen = cap(hcal, 19, [35,36,37,52,53], -412.9, 56.8, 324.6)
twenty = cap(solenoid_borders, 20, [54,55,56], -412.9, 348.3, 352.3)
twentyone = cap(solenoid_mid, 21, [57,58,59], -412.9, 364.9, 399.3)
twentytwo = cap(solenoid_borders, 22, [60,61,62], -412.9, 425, 429)
twentythree = barrel(solenoid_borders, 23, [54,56], 348.3, -412.9, 412.9)
twentyfour = barrel(solenoid_mid, 24, [57,59], 364.9, -412.9, 412.9)
twentyfive = barrel (solenoid_borders, 25, [60,62], 425, -412.9, 412.9)
twentysix = barrel(muon_detector, 26, [4,5], 446.1, -417.9, 417.9)
twentyseven = copy.deepcopy(nine)
twentyseven.id = 27
twentyseven.density = nozzle.density
twentyseven.next_ids = [9,10,11,12,13]
twentyeight = barrel(hcal, 28, [35,36,37,52,53], 324.6, -412.9, -235.4 )
twentynine = cap(hcal, 29, [63,64,65],-221, 174, 333)
thirty = barrel(hcal, 30, [63,64,65], 333, -221, 221)
thirtyone = cap(muon_detector, 31, [5,66], 417.9, 57.5, 446.1)
thirtytwo = barrel(solenoid_borders, 32, [54,55,56], 352.3, -412.9, 412.9)
thirtythree = barrel(solenoid_mid, 33, [57,58,59], 399.3, -412.9, 412.9)
thirtyfour = barrel(solenoid_borders, 34, [60,61,62], 429, -412.9, 412.9)

thirtyfive = cap(ecal, 35, [38,39,40,67], -250.9, 33.9, 170)
thirtysix = cap(space_between, 36, [29,68,69,70,71,72], -235.4, 170, 324.6)
thirtyseven = barrel(space_between, 37, [23,29],324.6, -412.9, -235.4 )

thirtyeight = cap(space_between, 38, [41,42,43,44,69, 70, 71, 72, 87], -230.7, 31, 170)
thirtynine = barrel(hcal, 39, [36], 170, -250.9, -235.4)
forty = barrel(space_between, 40, [23,29,69,71], 170, -235.4, -230.7)

fortyone = barrel(ecal, 41, [73,74], 150, -221, 221)
fortytwo = cap(ecal, 42, [75,76,77,78], 230.7, 31, 170)
fortythree = conic(nozzle, 43, [45,46,47,48,49,50,51], (31-2.2)/(230.7 - 6.5), 6.5, 230.7, 2.2, 'away')
fortyfour = barrel(beampipe, 44, [16,17,18], 2.2, -6.5, 6.5)
fortyfive = copy.deepcopy(fortythree)
fortyfive.id = 45
fortyfive.density = space_between.density
fortyfive.next_ids = [41,42]
seventyeight = conic(nozzle, 78, [46,47,49,50,51], (33.9 - 31)/(250.9 - 230.7), 230.7, 250.9, 31, 'away')
fortysix = copy.deepcopy(seventyeight)
fortysix.id = 46
fortysix.density = ecal.density
fortysix.next_ids = [77]
fortyseven = conic(hcal, 47, [79,80], (56.8 - 33.9)/(412.9 - 250.9), 250.9, 412.9, 33.9, 'away')
fortyeight = conic(space_between, 48, [31], (57.5 - 56.8)/( 417.9 - 412.9), 412.9, 417.9, 56.8, 'away')
fortynine = conic(muon_detector, 49, [5], (78.2 - 57.5)/(563.8 - 417.9), 417.9, 563.8, 57.5, 'away')
fifty = cap(outside, 50, [], 563.8, 2.2, 78.2)
fiftyone = barrel(beampipe, 51, [17,18], 2.2, 6.5, 563.8)

fiftytwo = barrel(ecal, 52, [38,67], 170, -250.9, -235.4)
fiftythree = copy.deepcopy(ten)
fiftythree.id = 53
fiftythree.density = nozzle.density
fiftythree.next_ids = [10,11,12,13]
fiftyfour = barrel(space_between, 54, [24,31], 352.3, -412.9, 412.9)
fiftyfive = barrel(space_between, 55, [23,28,29,30,31,81],  348.3, -412.9, 412.9)
fiftysix = cap(space_between, 56, [31], 412.9, 348.3, 352.3)
fiftyseven = barrel(space_between, 57, [25,31], 399.3, -412.9, 412.9)
fiftyeight = barrel(space_between, 58, [24,31,32], 364.9, -412.9, 412.9)
fiftynine = cap(space_between, 59, [31], 412.9, 364.9, 399.3)
sixty = barrel(space_between, 60, [26,31], 429, -412.9, 412.9)
sixtyone = barrel(space_between, 61, [25,31,33], 425, -412.9, 412.9)
sixtytwo = cap(space_between, 62, [31], 412.9,425, 429 )

sixtythree = barrel(space_between, 63, [23,31], 333, -221, 221)
sixtyfour = cap(space_between, 64, [23,31,72], 221, 174, 333)
sixtyfive = barrel(space_between, 65, [70,72,83], 174, -221, 221)
sixtysix = copy.deepcopy(fortynine)
sixtysix.id = 66
sixtysix.density = nozzle.density
sixtysix.next_ids = [49,50,51]
sixtyseven = copy.deepcopy(eleven)
sixtyseven.id = 67
sixtyseven.density = nozzle.density
sixtyseven.next_ids = [11,12,13]
sixtyeight = copy.deepcopy(forty)
sixtyeight.id = 68
sixtyeight.density = ecal.density
sixtyeight.next_ids = [35]
sixtynine = cap(ecal, 69, [73,74,84], -221, 150.0, 170.2)
seventy = barrel(ecal, 70, [73,74,84], 170.2, -221, 221)
seventyone = copy.deepcopy(sixtyfive)
seventyone.id = 71
seventyone.density = hcal.density
seventyone.next_ids = [63,64]
seventytwo = cap(hcal, 72, [79,80,85], 235.4, 170, 324.6)
seventythree = copy.deepcopy(seventy)
seventythree.id = 73
seventythree.density = space_between.density
seventythree.next_ids = [71,72]
seventyfour = cap(space_between, 74, [42,71,72], 221, 150, 170.2)
seventyfive = barrel(space_between, 75, [72], 170, 230.7, 235.4)
seventysix = barrel(hcal, 76, [79,80], 170, 235.4, 250.9)
seventyseven = cap(hcal, 77, [79,80,85], 250.9, 33.9, 170)
seventynine = cap(space_between, 79, [31,82], 412.9, 56.8, 324.6)
eighty = barrel(space_between, 80, [23,31], 324.6, 235.4, 412.9)
eightyone = copy.deepcopy(eighty)
eightyone.id = 81
eightyone.density = hcal.density
eightyone.next_ids = [79,80,85]
eightytwo = copy.deepcopy(fortyeight)
eightytwo.id = 82
eightytwo.density = nozzle.density
eightytwo.next_ids = [48,49,50,51]
eightythree = copy.deepcopy(seventyfive)
eightythree.id = 83
eightythree.density = ecal.density
eightythree.next_ids = [77]
eightyfour = copy.deepcopy(fortyone)
eightyfour.id = 84
eightyfour.density = space_between.density
eightyfour.next_ids = [41,42,43,44,87]
eightyfive = copy.deepcopy(fortyseven)
eightyfive.id = 85
eightyfive.density = nozzle.density
eightyfive.next_ids = [47,48,49,50,51]
eightysix = copy.deepcopy(seventysix)
eightysix.id = 86
eightysix.density = ecal.density
eightysix.next_ids = [77]
eightyseven = copy.deepcopy(twelve)
eightyseven.id = 87
eightyseven.density = nozzle.density
eightyseven.next_ids = [12,13]







OBJECTS = [zero, one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, 
           fifteen, sixteen, seventeen, eighteen, nineteen, twenty, twentyone, twentytwo, twentythree, twentyfour, 
           twentyfive, twentysix, twentyseven, twentyeight,twentynine,thirty, thirtyone, thirtytwo, thirtythree, thirtyfour, 
           thirtyfive, thirtysix, thirtyseven, thirtyeight, thirtynine, forty, fortyone, fortytwo, fortythree, fortyfour, fortyfive, 
           fortysix, fortyseven, fortyeight, fortynine, fifty, fiftyone, fiftytwo, fiftythree, fiftyfour, fiftyfive, fiftysix, fiftyseven,
           fiftyeight, fiftynine, sixty, sixtyone, sixtytwo, sixtythree, sixtyfour, sixtyfive, sixtysix, sixtyseven, sixtyeight,
           sixtynine, seventy, seventyone, seventytwo, seventythree, seventyfour, seventyfive, seventysix, seventyseven, seventyeight,
           seventynine, eighty, eightyone, eightytwo, eightythree, eightyfour, eightyfive, eightysix, eightyseven]
zbeginning = -563.8
rmax = 645
rbp = 2.2
zending = -1 * zbeginning
INITIALIZERS = [one, two, three] #should not contain minus_one or zero or minus_two
DECAYER = [minus_one, minus_two]
OUTSIDE = [minus_one, zero, four, five, fifty, eighteen]
iterations = 35
facedict = {'muon_detector': [1, 26, 31, 8, 49], 
            'solenoid': [20,21,22,23,24,25,32,33,34], 
            'hcal':[10,19,28,29,30,39,47,71,72,76,77,81], 
            'ecal':[11,35,41,42,46,52,68,69,70,83,86], 
            'nozzles':[-2,2,14,15,17,27,43,78,53,66,67,82,85,87]}
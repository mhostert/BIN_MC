
# approximate muon detector (v1)

from helpers import face, cap, barrel, conic, initializer, decayer
import useful_data as ud

#these densities will be number densities
block = face(19.3 / (939e6 *1.6e-19/(3e8)**2 * 10**3))
beampipe = face(0)
outside = face(0)

# caps have: id, next_ids, zpos, rbeg, rend

# barrels have: id, next_ids, rpos, zbeg, zend

# conic have: id, next_ids, tan_theta, zbeg, zend, rsmall, direction

# initializers have: id, next_ids, rsmall, rbig (their zpos is already defined in )

# decayers have: id, next_ids

minus_two = decayer(block, -2, [3,5,6])
minus_one = decayer(beampipe,-1, [3,4])

zero = cap(outside,0, [], 0,0,0) 
one = initializer(block, 1, [5,6,7],3, 645)
two = initializer(beampipe,2, [3,4], 0,3)

three = barrel(block, 3, [5,6], 3, -564, 564)
four = cap(outside, 4, [], 564, 0, 3)
five = barrel(outside, 5, [], 645, -564, 564)
six = cap(outside, 6, [], 564, 3, 645)
seven = barrel(beampipe, 7, [3,4], 3, -564, 564)


OBJECTS = [zero, one, two, three, four, five, six, seven]
zbeginning = -564
rmax = 645
rbp = 3
zending = -1 * zbeginning
INITIALIZERS = [one, two] #should not contain minus_one or zero
DECAYER = [minus_one, minus_two]
OUTSIDE = [zero,four, five, six]
iterations = 7
facedict = {}
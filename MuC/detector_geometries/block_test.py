# approximate muon detector (v1)
import MuC.detector_tools as det


# these densities will be number densities
block = det.Material(
    density=19.3 / (939e6 * 1.6e-19 / (3e8) ** 2 * 10**3), am=1, A=1, Z=1
)
beampipe = det.vacuum
outside = det.vacuum

# caps have: id, next_ids, zpos, rbeg, rend

# barrels have: id, next_ids, rpos, zbeg, zend

# conic have: id, next_ids, tan_theta, zbeg, zend, rsmall, direction

# initializers have: id, next_ids, rsmall, rbig (their zpos is already defined in )

# decayers have: id, next_ids

minus_two = det.MuonContainer(block, -2, [3, 5, 6], 6)
minus_one = det.MuonContainer(beampipe, -1, [3, 4], 4)

zero = det.BarrelCap(outside, 0, [], 0, 0, 0)
one = det.InitialFaces(block, 1, [5, 6, 7], 3, 645)
two = det.InitialFaces(beampipe, 2, [3, 4], 0, 3)

three = det.Barrel(block, 3, [5, 6], 3, -564, 564)
four = det.BarrelCap(outside, 4, [], 564, 0, 3)
five = det.Barrel(outside, 5, [], 645, -564, 564)
six = det.BarrelCap(outside, 6, [], 564, 3, 645)
seven = det.Barrel(beampipe, 7, [3, 4], 3, -564, 564)


OBJECTS = [zero, one, two, three, four, five, six, seven]
zbeginning = -564
rmax = 645
beam_pipe_radius = 3
zending = -1 * zbeginning
INITIALIZERS = [one, two]  # should not contain minus_one or zero
DECAYER = [minus_one, minus_two]
OUTSIDE = [zero, four, five, six]
iterations = 7
facedict = {}
name = "Det v0"

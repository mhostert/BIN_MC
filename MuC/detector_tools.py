import copy
import numpy as np

from numba import njit
from scipy.optimize import bisect

# from memory_profiler import profile

from DarkNews import Cfourvec as Cfv


AVOGADRO = 6.02214e23
MASS_NUCLEON = 939e6 * 1.6e-19 / (3e8) ** 2 * 10**3  # g
LIGHT_SPEED = 2.998e10  # cm/s


class face:
    """A detector component. Misnamed. Its density is the number of targets density."""

    def __init__(self, density):
        self.density = density  # density in which they're going


class cap:
    """These are the ones sitting on a z plane"""

    def __init__(
        self, parent, id, next_ids, zpos, rbeg, rend
    ):  # what to do for the end
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.zpos = zpos
        self.rbeg = rbeg
        self.rend = rend

    def check_intersection(
        self, position, momenta, mask
    ):  # position and momentum will be in (sample_size, 3) shape, 0,0,0 is center of detector
        return cap_check_i(self.zpos, self.rbeg, self.rend, position, momenta, mask)


class barrel:
    """These are cylindrical around a z line; usually both sides are included as separate faces"""

    def __init__(self, parent, id, next_ids, rpos, zbeg, zend):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.rpos = rpos
        self.zbeg = zbeg
        self.zend = zend

    def check_intersection(self, position, momenta, mask):
        indices = np.where(mask)[
            0
        ]  # array of indices of particles that we will consider
        a, b, c = barrel_get_pols(self.rpos, position[indices], momenta[indices])
        info = get_sols(
            a, b, c, self.zbeg, self.zend, position[mask], momenta[mask], mask
        )

        return info[0], info[1]


class conic:
    """These are cones centered on a z line"""

    def __init__(self, parent, id, next_ids, tan_theta, zbeg, zend, rsmall, direction):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.tan_theta = tan_theta  # positive
        self.rsmall = rsmall
        self.zbeg = zbeg  # big opening
        self.zend = zend  # smallest radius point
        self.direction = direction

        if self.direction == "towards":
            self.zcenter = zend + rsmall / tan_theta

        else:
            self.zcenter = zbeg - rsmall / tan_theta

    def check_intersection(self, position, momenta, mask):
        a, b, c = conic_get_pols(
            self.zcenter, self.tan_theta, position[mask], momenta[mask]
        )
        info = get_sols(
            a, b, c, self.zbeg, self.zend, position[mask], momenta[mask], mask
        )

        return info[0], info[1]


class initializer:
    """These are the first components a particle will encounter. Since they are all on the same plane (caps), I figured it's easier to compute once than to check at each time."""

    def __init__(self, parent, id, next_ids, rsmall, rbig):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.rsmall = rsmall
        self.rbig = rbig

    def check_in(self, r, mask):
        indices = np.where(mask)[0]
        new_mask = r[indices] < self.rbig

        return indices[np.where(new_mask)[0]]


class decayer:
    """These are the parts within which particles decaying in the detector would be. Supposed to only be beampipe, but in completely circular approximation might be nozzles."""

    def __init__(self, parent, id, next_ids, last):
        self.density = parent.density
        self.id = id
        self.next_ids = next_ids
        self.last = last

    def check_in(self, neighbor, zs, mask_decay, loc_mask):
        indices = np.where(mask_decay)[0]

        if neighbor.id == self.last:

            return indices[np.where(loc_mask)[0]], loc_mask

        local_mask = (zs > 0) & ((zs < neighbor.zend) & (loc_mask))

        return indices[np.where(local_mask)[0]], local_mask


@njit
def cap_check_i(zpos, rbeg, rend, position, momenta, mask):
    """Finds the multiplication parameter t by which a vector reaches a cap."""
    indices = np.where(mask)[0]  # array of indices of particles that we will consider
    delta_z = zpos - position[indices][:, 2]
    t = delta_z / momenta[indices][:, 2]
    ip = position[indices] + t[:, np.newaxis] * momenta[indices]  # interaction points
    r_values = np.sqrt(ip[:, 0] ** 2 + ip[:, 1] ** 2)
    mask_new = (
        (r_values < rend) & (r_values > rbeg) & (t > 0)
    )  # particles that have a correct ip
    kept_indices = indices[
        np.where(mask_new)[0]
    ]  # indexing the indices to get the number of the particles

    return kept_indices, ip[mask_new]


@njit
def barrel_get_pols(rpos, position, momenta):
    """Polynomial coefficients of a cylinder with a line."""
    a = (momenta[:, 0]) ** 2 + (momenta[:, 1]) ** 2
    b = 2 * (position[:, 0] * momenta[:, 0] + position[:, 1] * momenta[:, 1])
    c = position[:, 1] ** 2 + position[:, 0] ** 2 - rpos**2

    return a, b, c


@njit
def conic_get_pols(zcenter, tan_theta, position, momenta):
    """Polynomial coefficients of a cone with a line."""
    delta_z = zcenter - position[:, 2]
    a = -1 * momenta[:, 2] ** 2 * tan_theta**2 + momenta[:, 1] ** 2 + momenta[:, 0] ** 2
    b = (
        2 * position[:, 0] * momenta[:, 0]
        + 2 * position[:, 1] * momenta[:, 1]
        + tan_theta**2 * 2 * zcenter * momenta[:, 2]
        - tan_theta**2 * 2 * position[:, 2] * momenta[:, 2]
    )
    c = -1 * tan_theta**2 * delta_z**2 + position[:, 0] ** 2 + position[:, 1] ** 2

    return a, b, c


@njit
def get_roots(a, b, c):
    """Generates roots of polynomials represented by coefficients a,b,c."""
    disc = b**2 - 4 * a * c
    roots = np.empty((a.shape[0], 2))

    mask_1 = disc > 0
    mask_2 = disc == 0
    mask_3 = disc < 0

    div = 2 * a
    r = -1 * b / div
    r2 = np.sqrt(disc[mask_1]) / div[mask_1]

    roots[mask_1, 0] = r[mask_1] + r2
    roots[mask_1, 1] = r[mask_1] - r2
    roots[mask_2, 0] = r[mask_2]
    roots[mask_2, 1] = r[mask_2]
    roots[mask_3, :] = -1

    return roots


@njit
def get_sols(a, b, c, zbeg, zend, position, momenta, mask):
    """Based on polynomial solutions and characteristics of a component, finds which particles intersect it."""
    indices = np.where(mask)[0]
    roots = get_roots(a, b, c)
    ip_1 = position[:, 2] + roots[:, 0] * momenta[:, 2]  # only z
    ip_2 = position[:, 2] + roots[:, 1] * momenta[:, 2]  # only z

    new_mask_1 = (
        (ip_1 > zbeg) & (ip_1 < zend) & (np.round(roots[:, 0], decimals=12) > 0)
    )
    new_mask_2 = (
        (ip_2 > zbeg) & (ip_2 < zend) & (np.round(roots[:, 1], decimals=12) > 0)
    )
    t = roots.copy()  # need t to get the correct root
    any = new_mask_1 | new_mask_2
    doubles = new_mask_1 & new_mask_2
    t[new_mask_2, 0] = roots[new_mask_2, 1]

    if np.any(doubles):
        t[doubles, 0] = np.min(roots[doubles])

    ip = position[any] + t[any, 0][:, np.newaxis] * momenta[any]
    kept_indices = indices[
        np.where(any)[0]
    ]  # indexing the indices to get the correct particle numbers (ids)

    return kept_indices, ip


class cc:
    """Instance of a coordinate container (?). Holds all necessary information
    relating to the Monte Carlo muon decay simulation, and of the storage ring geometry
    """

    def __init__(self, C, w, sample_size, N_mu, pnumu, pnue, pos, name, Emu):
        self.name = name
        self.C = C
        nw = np.copy(w)
        self.weights = nw / np.sum(nw)
        self.sample_size = sample_size
        self.Nmu = N_mu
        self.p = pos
        self.pnumu = pnumu
        self.pnue = pnue
        self.vmu = 2.998e10 * np.sqrt(
            (1 - (105.6583755e-3 / Emu) ** 2)
        )  # speed of muons

    # @profile
    def straight_segment_at_detector(self, zlength, Lss, two=False):
        """Changes the coordinate axis, position, and momenta of particles to fit a storage ring geometry.
        Lss = -1 means that the cc object has already been transformed."""
        sim = copy.deepcopy(self)
        # print(sim.C)
        if Lss == -1:
            return sim

        elif Lss == 0:
            Lss = 0.01

        Lss = 100 * Lss
        sim.L = Lss
        K = Lss / (sim.C - Lss)
        a = 1e-6
        b = np.pi - 1e-6
        eta = bisect(to_opt, a, b, args=(K,))
        sim.Racc = Lss / 2 / np.sin(eta)

        factor = 0

        if two:
            factor = sim.C / 2

        zpos = (sim.p[:, 2] + factor) % sim.C

        if not Lss:
            on_straight_mask_left = [False] * sim.sample_size
            on_straight_mask_right = [False] * sim.sample_size

        else:
            on_straight_mask_left = zpos > (sim.C - Lss / 2)
            on_straight_mask_right = zpos < Lss / 2

        # for straight segment decays
        on_circ = ~((on_straight_mask_right) | (on_straight_mask_left))
        sim.p[on_straight_mask_right, 2] = zpos[on_straight_mask_right]
        sim.p[on_straight_mask_left, 2] = zpos[on_straight_mask_left] - sim.C

        # lower dim quantities (for circ)
        phis = (zpos[on_circ] - Lss / 2) / sim.Racc + eta
        sim.pnumu[on_circ, :] = Cfv.rotationx(sim.pnumu[on_circ], phis)
        sim.pnue[on_circ, :] = Cfv.rotationx(sim.pnue[on_circ], phis)
        sim.p[on_circ, 2] = (sim.Racc + sim.p[on_circ, 1]) * np.sin(phis)
        sim.p[on_circ, 1] = (sim.Racc + sim.p[on_circ, 1]) * np.cos(
            phis
        ) - Lss / 2 / np.tan(eta)

        # for mltd
        mltd = np.empty(sim.sample_size)
        mltd[zpos < sim.C / 2] = -1 * zpos[zpos < sim.C / 2]
        mltd[zpos > sim.C / 2] = sim.C - zpos[zpos > sim.C / 2]

        # free memory
        mask_acc = (sim.p[:, 2] < zlength) & (sim.p[:, 1] > -1 * Lss / 2 / np.tan(eta))

        # NOTE: there are two masks, but one is never used? Why?
        mask_acc2 = (sim.p[:, 2] > -1 * zlength) & (
            sim.p[:, 1] > -1 * Lss / 2 / np.tan(eta)
        )

        sim.p = sim.p[mask_acc]
        sim.pnumu = sim.pnumu[mask_acc]
        sim.pnue = sim.pnue[mask_acc]
        sim.weights = sim.weights[mask_acc]
        sim.sample_size = np.sum(mask_acc)
        sim.times = mltd[mask_acc] / self.vmu[mask_acc]

        if two:
            sim2, _ = self.straight_segment_at_detector(
                zlength, Lss=Lss / 100, two=False
            )
            return sim, sim2

        return sim, None

    def clear_mem(self):
        """Freeing memory at the end of a sim, if necessary"""
        deletables = ["Racc", "p", "pnumu", "pnue", "weights", "sample_size", "Nmu"]
        for att in deletables:

            if hasattr(self, att):
                delattr(self, att)


class material:
    """To store material information for detector component densities and cross sections."""

    def __init__(self):
        self.N = 0
        self.density = 0
        self.e = 0


class subs(material):
    """Pure substances; periodic elements."""

    def __init__(self, density, am, A, Z):
        super().__init__()
        self.density = density
        self.am = am
        self.Z = Z
        self.A = A
        nq = AVOGADRO * self.density / self.am
        self.N = nq * self.A
        self.e = nq * self.Z


class comp(material):
    """Compositions of periodic elements."""

    def __init__(self, table):
        super().__init__()

        for row in table:
            self.density += row[0].density * row[1]
            self.N += (
                row[0].N * row[1]
            )  # row[0] is N of an element; row[1] is the percentage of it that occupies the total material
            self.e += row[0].e * row[1]


class unif(material):
    """Alloys/compositions of uniform densities."""

    def __init__(self, density):
        super().__init__()
        self.density = density
        self.N = self.density / MASS_NUCLEON
        self.e = self.N / 2


def to_opt(x, K):
    """Optimizing provides root for transformed radius, R', of storage ring for nonzero Lss."""
    return np.sin(x) - K * (np.pi - x)


# Pre-defined substances

# density in g/cm**3; atomic mass in g/mol
Si = subs(2.329, 28.0855, 28, 14)
WSi2 = subs(9.3, 240.01, 240, 102)
Fe = subs(7.874, 55.845, 56, 26)
Al = subs(2.7, 26.981539, 27, 13)
W = subs(19.3, 183.84, 184, 74)
Cu = subs(8.96, 63.546, 64, 29)
PS = subs(1.05, 104.1, 104, 56)

# from CLICdet paper
hcal_CLICdet = comp(
    [[Fe, 20 / 26.5], [Al, 0.7 / 26.5], [Cu, 0.1 / 26.5], [PS, 3 / 26.5]]
)
ecal_CLICdet = comp([[W, 1.9 / 5.05], [Cu, 2.3 / 5.05], [Si, 0.5 / 5.05]])

# from online
EARTH_DENSITY = unif(5.51).N
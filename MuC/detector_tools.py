import numpy as np
from numba import njit

from DarkNews import const

# from memory_profiler import profile


class Component:
    """A detector component. Its density is the number of targets density."""

    def __init__(self, density):
        self.density = density  # density in which they're going


class BarrelCap:
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


class Barrel:
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
        a, b, c = barrel_get_polynomial(self.rpos, position[indices], momenta[indices])
        info = get_sols(
            a, b, c, self.zbeg, self.zend, position[mask], momenta[mask], mask
        )

        return info[0], info[1]


class Conic:
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
        a, b, c = conic_get_polynomial(
            self.zcenter, self.tan_theta, position[mask], momenta[mask]
        )
        info = get_sols(
            a, b, c, self.zbeg, self.zend, position[mask], momenta[mask], mask
        )

        return info[0], info[1]


class InitialFaces:
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


class MuonContainer:
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
def barrel_get_polynomial(rpos, position, momenta):
    """Polynomial coefficients of a cylinder with a line."""
    a = (momenta[:, 0]) ** 2 + (momenta[:, 1]) ** 2
    b = 2 * (position[:, 0] * momenta[:, 0] + position[:, 1] * momenta[:, 1])
    c = position[:, 1] ** 2 + position[:, 0] ** 2 - rpos**2

    return a, b, c


@njit
def conic_get_polynomial(zcenter, tan_theta, position, momenta):
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


class Material:
    """Pure substances; periodic elements."""

    def __init__(self, density, am, A, Z):
        self.density = density
        self.am = am
        self.Z = Z
        self.A = A
        nq = const.NAvo * self.density / self.am
        self.N = nq * self.A
        self.e = nq * self.Z


class CompositMaterial:
    def __init__(self, table):
        """Compositions of materials.
        fraction is the percentage of it that occupies the total material
        """
        self.density = 0
        self.N = 0
        self.e = 0
        for material, fraction in table:
            self.density += material.density * fraction
            self.N += material.N * fraction
            self.e += material.e * fraction


# class unif(Material):
#     """Alloys/compositions of uniform densities."""

#     def __init__(self, density):
#         super().__init__()
#         self.density = density
#         self.N = self.density / (const.m_avg / const.g_to_GeV)
#         self.e = self.N / 2


# Pre-defined substances

# density in g/cm**3; atomic mass in g/mol
Si = Material(2.329, 28.0855, 28, 14)
WSi2 = Material(9.3, 240.01, 240, 102)
Fe = Material(7.874, 55.845, 56, 26)
Al = Material(2.7, 26.981539, 27, 13)
W = Material(19.3, 183.84, 184, 74)
Cu = Material(8.96, 63.546, 64, 29)
PS = Material(1.05, 104.1, 104, 56)

# from CLICdet paper
hcal_CLICdet = CompositMaterial(
    [[Fe, 20 / 26.5], [Al, 0.7 / 26.5], [Cu, 0.1 / 26.5], [PS, 3 / 26.5]]
)
ecal_CLICdet = CompositMaterial([[W, 1.9 / 5.05], [Cu, 2.3 / 5.05], [Si, 0.5 / 5.05]])

OinAir = Material(1.225e-3, 15.9994, 16, 8)  # air density in g/cm**3
NinAir = Material(1.225e-3, 14.0067, 14, 7)
ArinAir = Material(1.225e-3, 39.95, 40, 18)
Air = CompositMaterial(
    [
        [OinAir, 0.2095],
        [NinAir, 0.7812],
        [ArinAir, 0.0093],
    ]
)

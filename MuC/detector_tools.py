import numpy as np
from numba import njit
<<<<<<< HEAD

from DarkNews import const
=======
from scipy.optimize import bisect
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929

# from memory_profiler import profile


<<<<<<< HEAD
class Component:
    """A detector component. Its density is the number of targets density."""
=======

AVOGADRO = 6.02214e23
MASS_NUCLEON = 939e6 * 1.6e-19 / (3e8) ** 2 * 10**3  # g
LIGHT_SPEED = 2.998e10  # cm/s
ry_max = 1434.4  # m
total_arc_length = 9851.838522846553


class face:
    """A detector component. Misnamed. Its density is the number of targets density."""
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929

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


<<<<<<< HEAD
class Material:
=======
class cc:
    """Instance of a coordinate container (?). Holds all necessary information
    relating to the Monte Carlo muon decay simulation, and of the storage ring geometry
    """

    def __init__(self, C, w, sample_size, Nmu_per_bunch, pnumu, pnue, pos, name, Emu):
        self.name = name
        self.C = C
        nw = np.copy(w)
        self.weights = nw / np.sum(nw)
        self.sample_size = sample_size
        self.Nmu_per_bunch = Nmu_per_bunch
        self.p = pos
        self.pnumu = pnumu
        self.pnue = pnue
        self.vmu = 2.998e10 * np.sqrt(
            (1 - (105.6583755e-3 / Emu) ** 2)
        )  # speed of muons

    # @profile
    def straight_segment_at_detector(self, zlength, Lss, two=False):
        """Changes the coordinate axis, position, and momenta of particles to fit a storage ring geometry.
        zlength is the length of the detector on one side.
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

    def parametrized_curve(self, zlength, Lss, two=False, smoother=25):
        # already have self.name, self.C (need to change), self.weights, self.sample_size, self.Nmu, self.p (position), self.pnumu (4M), self.pnue, self.vmu
        self.C = total_arc_length * 100

        sim = copy.deepcopy(self)

        if Lss == -1:
            return sim
        elif not ((Lss <= 400) & (Lss >= 75)):
            raise ValueError(
                "The only possible striaght segment values for the parametrized curve are 75-400 m."
            )
        Lss = 100 * Lss  # cm
        sim.L = Lss
        factor = 0

        if two:
            factor = sim.C / 2

        t_values = (sim.p[:, 2] + factor) % sim.C  # not translated yet

        # reduces computation; nothing to do with straight segment yet
        mask = (t_values < total_arc_length * 100 / 4) | (t_values > sim.C -1 * zlength * 2)

        # straight segment stuff is hidden within parametric_position_by_t_array; Lss is in cm already
        sim.p[mask, :] = parametric_position_by_t_array(t_values[mask], Lss, smoother)

        phis = get_theta_for_p_rotation(t_values[mask], Lss, smoother)

        sim.pnumu[mask, :] = Cfv.rotationx(sim.pnumu[mask], -1 * phis)
        sim.pnue[mask, :] = Cfv.rotationx(sim.pnue[mask], -1 * phis)

        # for mltd
        mltd = np.empty(sim.sample_size)
        mltd[t_values < sim.C / 2] = t_values[t_values < sim.C / 2]
        mltd[t_values > sim.C / 2] = -1*sim.C + t_values[t_values > sim.C / 2]

        mask_acc = mask
        sim.p = sim.p[mask_acc]
        sim.pnumu = sim.pnumu[mask_acc]
        sim.pnue = sim.pnue[mask_acc]
        sim.weights = sim.weights[mask_acc]
        sim.sample_size = np.sum(mask_acc)
        sim.times = mltd[mask_acc] / self.vmu[mask_acc]

        if two:
            sim2, _ = self.parametrized_curve(zlength, Lss=Lss / 100, two=False)
            return sim, sim2

        return sim, None

    def clear_mem(self):
        """Freeing memory at the end of a sim, if necessary"""
        deletables = [
            "Racc",
            "p",
            "pnumu",
            "pnue",
            "weights",
            "sample_size",
            "Nmu_per_bunch",
        ]
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
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929
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

<<<<<<< HEAD
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
=======
# from online
EARTH_DENSITY = unif(5.51).N


# Precompute func_radius for all theta values
def func_radius(theta):
    """This is the regular parametrization; straight segment comes after. DO NOT CHANGE unless altering the actual curve parametrization."""
    radius = 1630
    matching_angular_distance = 25 / 180 * np.pi
    tan_theta = np.abs(np.tan(theta))
    condition = tan_theta > np.tan(matching_angular_distance)
    arc_tan_value = np.where(condition, np.arctan(tan_theta), matching_angular_distance)
    result = (
        (arc_tan_value - matching_angular_distance)
        / (np.pi / 2 - matching_angular_distance)
        * (np.pi / 2)
    )
    return radius + (-radius * 1.2 / 10) * np.sin(result) ** 2


# Derivatives of x(θ) and y(θ) for arc length calculation
def dx_dtheta(theta, func_radius_vals):
    return -func_radius_vals * np.sin(theta) + np.cos(theta) * np.gradient(
        func_radius_vals, theta
    )


def dy_dtheta(theta, func_radius_vals):
    return func_radius_vals * np.cos(theta) + np.sin(theta) * np.gradient(
        func_radius_vals, theta
    )


theta = np.linspace(0, 2 * np.pi, 5000)  # Increase resolution for smooth interpolation
func_radius_vals = func_radius(theta)
# Calculate ds, arc length as a function of θ
ds = np.sqrt(
    dx_dtheta(theta, func_radius_vals) ** 2 + dy_dtheta(theta, func_radius_vals) ** 2
)
arc_length = cumtrapz(ds, theta, initial=0)

# Interpolate to get θ from arc length
theta_of_s = interp1d(
    arc_length, theta, bounds_error=False, fill_value="extrapolate"
)  # this input is in meters!


# Optimized function that takes a 1D array of t values and returns an array of coordinates
def parametric_position_by_t_array(t_values, Lss, smoother=25):

    t_values += total_arc_length * 100 / 4  # to start it at 0,0,0
    # print(arc_length[-1]) for total arc length
    # Ensure t_values are within the correct range [0, total_arc_length]; they are provided as cm though! Internally, deal with them as m.

    t_values_clipped = t_values / 100 % total_arc_length  # this is in meters

    # Get corresponding θ values for the given t_values
    theta_t = theta_of_s(t_values_clipped)  # this input is in meters!

    # Compute x_t and y_t for each θ
    x_t = func_radius(theta_t) * np.cos(theta_t)
    y_t = func_radius(theta_t) * np.sin(theta_t)

    # let's make the straight section here (flattening out the curve; overdense by about 5 cm for Lss = 200 m; more study about the actual effect of the contraction at larger LSS (400m) needs to be done.):

    # Note that these thetas take t = 0 to be pi/2.
    a1, a2, Rr, b1, b2 = sss(Lss, smoother)

    mask_ss = (theta_t < a2) & (theta_t > a1)

    # Translation by ry_max (assumed to be constant)
    result = np.column_stack([np.zeros_like(x_t), 100 * (y_t - ry_max), 100 * x_t])

    result[mask_ss, 1] = 0
    result[~mask_ss, 1] -= Rr * 100

    # smoothening mask
    mask_sm = (theta_t > b1) & (theta_t < b2)
    cc = smooth(theta_t[mask_sm], Lss, Rr, b1, b2, smoother)
    result[mask_sm, 1] = cc[1]
    result[mask_sm, 2] = cc[0]

    return result


def get_theta_for_p_rotation(t_values, Lss, smoother=25, dt=1):
    # note that the t parameter goes counterclockwise; our muons go clockwise. So, our momentum rotation needs to be treated carefully.
    xf = parametric_position_by_t_array(t_values - dt / 2, Lss, smoother)
    xi = parametric_position_by_t_array(t_values + dt / 2, Lss, smoother)
    phi = np.arctan2(
        (xf[:, 1] - xi[:, 1]), (xf[:, 2] - xi[:, 2])
    )  # these are the rotation phis (rotation clockwise about x) for momentum.
    return phi


def smooth(thetas, Lss, Rr, b1, b2, smoother=25):
    """Have to generalize this. How?..."""
    yf = (func_radius(b2) * np.sin(b2) - ry_max) * 100 - Rr * 100
    L = Lss / 2
    a = yf / 4 / (smoother * 100) ** 2
    b = 2 * a * (-1 * smoother * 100 + L)
    c = a * (L - smoother * 100) ** 2
    x0 = func_radius(b1) * np.cos(b1) * 100
    xf = func_radius(b2) * np.cos(b2) * 100
    interpx = np.linspace(x0, xf, 100)
    interpy = a * interpx**2 + b * interpx + c
    smoothener = interp1d(
        interpx, interpy, kind="linear", bounds_error=False, fill_value="extrapolate"
    )
    xs = (thetas - b1) / (b2 - b1) * (xf - x0) + x0
    ys = smoothener(xs)
    return xs, ys


def sss(Lss, smoother=25):
    """Straight segment shenanigans; input Lss is in cm"""

    a1 = theta_of_s(-Lss / 200 + total_arc_length / 4)
    a2 = theta_of_s(Lss / 200 + total_arc_length / 4)
    b1 = theta_of_s(Lss / 200 - smoother + total_arc_length / 4)
    b2 = theta_of_s(Lss / 200 + smoother + total_arc_length / 4)
    Rr = func_radius(a2) * np.sin(a2) - ry_max

    return a1, a2, Rr, b1, b2
>>>>>>> afbcf968c3e61d22c8bb0afe1be2991cd23da929

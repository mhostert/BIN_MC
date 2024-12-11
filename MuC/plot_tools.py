import numpy as np
from scipy.stats import chi2

from math import log10, floor, erf

import matplotlib
import colorsys
import matplotlib.colors as mc
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.tri as tri
from scipy.spatial.distance import pdist, squareform
from matplotlib import colors as mpl_colors
from matplotlib.collections import PatchCollection

import scipy

###########################
# Matheus
fsize = 11
fsize_annotate = 10

std_figsize = (1.2 * 3.7, 1.3 * 2.3617)
std_axes_form = [0.18, 0.16, 0.79, 0.76]

rcparams = {
    "axes.labelsize": fsize,
    "xtick.labelsize": fsize,
    "ytick.labelsize": fsize,
    "figure.figsize": std_figsize,
    "legend.frameon": False,
    "legend.loc": "best",
}
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}\usepackage{amssymb}"
rc("text", usetex=True)
rc("font", **{"family": "serif", "serif": ["Computer Modern Roman"]})
matplotlib.rcParams["hatch.linewidth"] = 0.3

rcParams.update(rcparams)

cblind_safe_wheel = [
    "#3f90da",
    "#ffa90e",
    "#bd1f01",
    "#94a4a2",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
]


##########################
#
def get_CL_from_sigma(sigma):
    return erf(sigma / np.sqrt(2))


def get_chi2vals_w_nsigmas(n_sigmas, ndof):
    return [chi2.ppf(get_CL_from_sigma(i), ndof) for i in range(n_sigmas + 1)]


def get_chi2vals_w_sigma(sigma, ndof):
    return chi2.ppf(get_CL_from_sigma(sigma), ndof)


def get_chi2vals_w_CL(CLs, ndof):
    return [chi2.ppf(cl, ndof) for cl in CLs]


def std_fig(ax_form=std_axes_form, figsize=std_figsize, rasterized=False):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(ax_form, rasterized=rasterized)
    ax.patch.set_alpha(0.0)
    return fig, ax


def double_axes_fig(
    height=0.5,
    gap=0.1,
    axis_base=[0.14, 0.1, 0.80, 0.18],
    figsize=std_figsize,
    split_y=False,
    split_x=False,
    rasterized=False,
):
    fig = plt.figure(figsize=figsize)

    if split_y and not split_x:
        axis_base = [0.14, 0.1, 0.80, 0.4 - gap / 2]
        axis_appended = [0.14, 0.5 + gap / 2, 0.80, 0.4 - gap / 2]

    elif not split_y and split_x:
        axis_appended = [0.14, 0.1, 0.4 - gap / 2, 0.8]
        axis_base = [0.14 + 0.4 + gap / 2, 0.1, 0.4 - gap / 2, 0.8]

    else:
        axis_base[-1] = height
        axis_appended = axis_base + np.array(
            [0, height + gap, 0, 1 - 2 * height - gap - axis_base[1] - 0.07]
        )

    ax1 = fig.add_axes(axis_appended, rasterized=rasterized)
    ax2 = fig.add_axes(axis_base, rasterized=rasterized)
    ax1.patch.set_alpha(0.0)
    ax2.patch.set_alpha(0.0)

    return fig, [ax1, ax2]


def data_plot(ax, X, Y, xerr, yerr, zorder=2, label="data", **kwargs):
    return ax.errorbar(
        X,
        Y,
        yerr=yerr,
        xerr=xerr,
        marker="o",
        markeredgewidth=0.75,
        capsize=1,
        markerfacecolor="black",
        markeredgecolor="black",
        ms=1.75,
        lw=0.0,
        elinewidth=0.75,
        color="black",
        label=label,
        zorder=zorder,
        **kwargs,
    )


def step_plot(
    ax, x, y, lw=1, color="red", label="signal", where="post", dashes=(3, 0), zorder=3
):
    return ax.step(
        np.append(x, np.max(x) + x[-1]),
        np.append(y, 0.0),
        where=where,
        lw=lw,
        dashes=dashes,
        color=color,
        label=label,
        zorder=zorder,
    )


def plot_MB_vertical_region(ax, color="dodgerblue", label=r"MiniBooNE $1 \sigma$"):
    ##########
    # MINIBOONE 2018
    matplotlib.rcParams["hatch.linewidth"] = 0.7
    y = [0, 1e10]
    NEVENTS = 381.2
    ERROR = 85.2
    xleft = (NEVENTS - ERROR) / NEVENTS
    xright = (NEVENTS + ERROR) / NEVENTS
    ax.fill_betweenx(
        y,
        [xleft, xleft],
        [xright, xright],
        zorder=3,
        ec=color,
        fc="None",
        hatch="\\\\\\\\\\",
        lw=0,
        label=label,
    )

    ax.vlines(1, 0, 1e10, zorder=3, lw=1, color=color)
    ax.vlines(xleft, 0, 1e10, zorder=3, lw=0.5, color=color)
    ax.vlines(xright, 0, 1e10, zorder=3, lw=0.5, color=color)


# Kevin align
def flushalign(ax):
    ic = 0
    for l in ax.get_yticklabels():
        if ic == 0:
            l.set_va("bottom")
        elif ic == len(ax.get_yticklabels()) - 1:
            l.set_va("top")
        ic += 1

    ic = 0
    for l in ax.get_xticklabels():
        if ic == 0:
            l.set_ha("left")
        elif ic == len(ax.get_xticklabels()) - 1:
            l.set_ha("right")
        ic += 1


# Function to find the path that connects points in order of closest proximity
def nearest_neighbor_path(points):
    # Compute the pairwise distance between points
    dist_matrix = squareform(pdist(points))

    # Set diagonal to a large number to avoid self-loop
    np.fill_diagonal(dist_matrix, np.inf)

    # Start from the first point
    current_point = 0
    path = [current_point]

    # Find the nearest neighbor of each point
    while len(path) < len(points):
        # Find the nearest point that is not already in the path
        nearest = np.argmin(dist_matrix[current_point])
        # Add the nearest point to the path
        path.append(nearest)
        # Update the current point
        current_point = nearest
        # Mark the visited point so it's not revisited
        dist_matrix[:, current_point] = np.inf

    # Return the ordered path indices and the corresponding points
    ordered_points = points[path]
    return ordered_points


def get_ordered_closed_region(points, logx=False, logy=False):
    xraw, yraw = points

    # check for nans
    if np.isnan(points).sum() > 0:
        raise ValueError("NaN's were found in input data. Cannot order the contour.")

    # check for repeated x-entries -- remove them
    # x, mask_diff = np.unique(x, return_index=True)
    # y = y[mask_diff]

    if logy:
        if (yraw == 0).any():
            raise ValueError("y values cannot contain any zeros in log mode.")
        yraw = np.log10(yraw)
    if logx:
        if (xraw == 0).any():
            raise ValueError("x values cannot contain any zeros in log mode.")
        xraw = np.log10(xraw)

    # Transform to unit square space:
    xmin, xmax = np.min(xraw), np.max(xraw)
    ymin, ymax = np.min(yraw), np.max(yraw)

    x = (xraw - xmin) / (xmax - xmin)
    y = (yraw - ymin) / (ymax - ymin)

    points = np.array([x, y]).T
    # points_s     = (points - points.mean(0))
    # angles       = np.angle((points_s[:,0] + 1j*points_s[:,1]))
    # points_sort  = points_s[angles.argsort()]
    # points_sort += points.mean(0)

    # if np.isnan(points_sort).sum()>0:
    #     raise ValueError("NaN's were found in sorted points. Cannot order the contour.")
    # # print(points.mean(0))
    # # return points_sort
    # tck, u = splprep(points_sort.T, u=None, s=0.0, per=0, k=1)
    # # u_new = np.linspace(u.min(), u.max(), len(points[:,0]))
    # x_new, y_new = splev(u, tck, der=0)
    # # x_new, y_new = splev(u_new, tck, der=0)
    dist_matrix = squareform(pdist(points))

    # Set diagonal to a large number to avoid self-loop
    np.fill_diagonal(dist_matrix, np.inf)

    # Start from the first point
    current_point = 0
    path = [current_point]

    # Find the nearest neighbor of each point
    while len(path) < len(points):
        # Find the nearest point that is not already in the path
        nearest = np.argmin(dist_matrix[current_point])
        # Add the nearest point to the path
        path.append(nearest)
        # Update the current point
        current_point = nearest
        # Mark the visited point so it's not revisited
        dist_matrix[:, current_point] = np.inf

    # Return the ordered path indices and the corresponding points
    x_new, y_new = points[path].T

    x_new = x_new * (xmax - xmin) + xmin
    y_new = y_new * (ymax - ymin) + ymin

    if logx:
        x_new = 10 ** (x_new)
    if logy:
        y_new = 10 ** (y_new)
    return x_new, y_new


def interp_grid(
    x,
    y,
    z,
    fine_gridx=False,
    fine_gridy=False,
    logx=False,
    logy=False,
    method="interpolate",
    smear_stddev=False,
):
    # default
    if not fine_gridx:
        fine_gridx = 100
    if not fine_gridy:
        fine_gridy = 100

    # log scale x
    if logx:
        xi = np.geomspace(np.min(x), np.max(x), fine_gridx)
    else:
        xi = np.linspace(np.min(x), np.max(x), fine_gridx)

    # log scale y
    if logy:
        yi = np.geomspace(np.min(y), np.max(y), fine_gridy)
    else:
        yi = np.linspace(np.min(y), np.max(y), fine_gridy)

    Xi, Yi = np.meshgrid(xi, yi)
    # if logy:
    #     Yi = 10**(-Yi)

    # triangulation
    if method == "triangulation":
        triang = tri.Triangulation(x, y)
        interpolator = tri.LinearTriInterpolator(triang, z)
        Zi = interpolator(Xi, Yi)

    elif method == "interpolate":
        Zi = scipy.interpolate.griddata(
            (x, y), z, (xi[None, :], yi[:, None]), method="linear", rescale=True
        )
    else:
        print(f"Method {method} not implemented.")

    # gaussian smear -- not recommended
    if smear_stddev:
        Zi = scipy.ndimage.filters.gaussian_filter(
            Zi, smear_stddev, mode="nearest", order=0, cval=0
        )

    return Xi, Yi, Zi


def round_sig(x, sig):
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def sci_notation(
    num,
    sig_digits=1,
    precision=None,
    exponent=None,
    notex=False,
    optional_sci=False,
):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if num != 0:
        if exponent is None:
            exponent = int(np.floor(np.log10(abs(num))))
        coeff = round(num / float(10**exponent), sig_digits)
        if coeff == 10:
            coeff = 1
            exponent += 1
        if precision is None:
            precision = sig_digits

        if optional_sci and np.abs(exponent) < optional_sci:
            string = rf"{round_sig(num, precision)}"
        else:
            string = r"{0:.{2}f}\times 10^{{{1:d}}}".format(coeff, exponent, precision)

        if notex:
            return string
        else:
            return f"${string}$"

    else:
        return r"0"


# https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


###########################
def get_cmap_colors(name, ncolors, cmin=0, cmax=1, reverse=False):
    try:
        cmap = plt.get_cmap(name)
    except ValueError:
        cmap = build_cmap(name, reverse=reverse)
    return cmap(np.linspace(cmin, cmax, ncolors, endpoint=True))


def build_cmap(color, reverse=False):
    cvals = [0, 1]
    colors = [color, "white"]
    if reverse:
        colors = colors[::-1]

    norm = plt.Normalize(min(cvals), max(cvals))
    tuples = list(zip(map(norm, cvals), colors))
    return mpl_colors.LinearSegmentedColormap.from_list("", tuples)


# define an object that will be used by the legend
class MulticolorPatch(object):
    def __init__(self, colors):
        self.colors = colors


# define a handler for the MulticolorPatch object
class MulticolorPatchHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        width, height = handlebox.width, handlebox.height
        patches = []
        for i, c in enumerate(orig_handle.colors):
            patches.append(
                plt.Rectangle(
                    [
                        width / len(orig_handle.colors) * i - handlebox.xdescent,
                        -handlebox.ydescent,
                    ],
                    width / len(orig_handle.colors),
                    height,
                    facecolor=c,
                    edgecolor="none",
                )
            )

        patch = PatchCollection(patches, match_original=True)

        handlebox.add_artist(patch)
        return patch


def plot_det(geom, ax, orientation="z-y", xl=True, yl=True):
    """Plots the detector geometry behind a sim plot."""

    if geom == "det_v1":
        T1 = [[-231, 231, 231, 28, -28, -231, -231], [150, 150, 24, 3, 3, 24, 150]]

        ECAL1 = [[-231, -231, 231, 231, -231], [150, 170, 170, 150, 150]]
        ECAL2 = [[231, 231, 251, 251, 231], [24, 170, 170, 26, 24]]
        ECAL3 = [[-1 * i for i in ECAL2[0]], ECAL2[1]]

        HCAL1 = [[-251, -251, 251, 251, -251], [170, 348, 348, 170, 170]]
        HCAL2 = [[251, 251, 418, 418, 251, 251, 251], [170, 348, 348, 43, 26, 170, 170]]
        HCAL3 = [[-1 * i for i in HCAL2[0]], HCAL2[1]]

        SOLENOID = [[-418, -418, 418, 418, -418], [348, 446, 446, 348, 348]]

        MD1 = [[-418, -418, 418, 418, -418], [446, 645, 645, 446, 446]]
        MD2 = [[418, 418, 564, 564, 418], [43, 645, 645, 58, 43]]
        MD3 = [[-1 * i for i in MD2[0]], MD2[1]]

        CONE1 = [[28, 564, 564, 28], [3, 58, 3, 3]]
        CONE2 = [[-1 * i for i in CONE1[0]], CONE1[1]]

        BL = [[-564, -564, 564, 564, -564], [-3, 3, 3, -3, -3]]

        components = [
            T1,
            ECAL1,
            ECAL2,
            ECAL3,
            HCAL1,
            HCAL2,
            HCAL3,
            SOLENOID,
            MD1,
            MD2,
            MD3,
            CONE1,
            CONE2,
            BL,
        ]
        cols = (
            ["lightgrey"] * 1
            + ["dimgrey"] * 3
            + ["grey"] * 3
            + ["darkgrey"]
            + 3 * ["grey"]
            + ["black"] * 2
            + ["white"]
        )

        for i, component in enumerate(components):
            ax.plot(component[0], component[1], color=cols[i], lw=0.5)
            ax.fill_between(component[0], component[1], color=cols[i], alpha=0.7)
            new_y = [-1 * k for k in component[1]]
            ax.plot(component[0], new_y, color=cols[i], lw=0.5)
            ax.fill_between(component[0], new_y, color=cols[i], alpha=0.7)

        ax.set_xlabel("z-coordinate (cm)")
        ax.set_ylabel("r-coordinate (cm)")

    elif (geom == "det_v2") | (geom == "zero_density_test"):

        if orientation == "x-y" or orientation == "y-x":
            MDET = 645
            SPS1 = 446.1
            SOL_1 = 429
            SPS2 = 425
            SOL_2 = 399.3
            SPS3 = 364.9
            SOL_3 = 352.3
            SPS4 = 348.3
            HCAL = 333
            SPS5 = 174
            ECAL = 170.2
            SPS6 = 150
            CONE = 31
            BL = 2.2
            components = [
                MDET,
                SPS1,
                SOL_1,
                SPS2,
                SOL_2,
                SPS3,
                SOL_3,
                SPS4,
                HCAL,
                SPS5,
                ECAL,
                SPS6,
                CONE,
                BL,
            ]
            cols = [
                "grey",
                "white",
                "gray",
                "white",
                "lightgrey",
                "white",
                "gray",
                "white",
                "darkgrey",
                "white",
                "dimgrey",
                "white",
                "black",
                "white",
            ]

            for i, component in enumerate(components):
                circle = plt.Circle(
                    (0, 0),
                    component,  # NOTE: det??
                    zorder=i / len(components),
                    alpha=1,
                    edgecolor=cols[i],
                    facecolor=cols[i],
                )
                ax.add_artist(circle)

            ax.scatter(MDET, MDET, color=cols[3])
            ax.scatter(MDET, MDET, color=cols[3])
            ax.scatter(MDET, MDET, color=cols[3])
            ax.scatter(MDET, MDET, color=cols[3])

            if orientation == "y-x":
                ax.set_ylim(-1 * 645 * 10 / 12, 645 * 10 / 12)
                ax.set_xlim(0, 645)
            elif orientation == "x-y":
                ax.set_xlim(-1 * 645 * 10 / 12, 645 * 10 / 12)
                ax.set_ylim(0, 645)

            if xl:
                ax.set_xlabel("x-coordinate (cm)")

            if yl:
                ax.set_ylabel("y-coordinate (cm)")

        else:
            ECAL1 = [[-221, -221, 221, 221, -221], [150, 170.2, 170.2, 150, 150]]
            ECAL2 = [[230.7, 230.7, 250.9, 250.9, 230.7], [31, 170, 170, 33.9, 31]]
            ECAL3 = [[-1 * i for i in ECAL2[0]], ECAL2[1]]

            HCAL1 = [[-221, -221, 221, 221, -221], [174, 333, 333, 174, 174]]
            HCAL2 = [
                [235.4, 235.4, 412.9, 412.9, 250.9, 250.9, 235.4],
                [170, 324.6, 324.6, 56.8, 33.9, 170, 170],
            ]
            HCAL3 = [[-1 * i for i in HCAL2[0]], HCAL2[1]]

            SOLENOID = [
                [-412.9, -412.9, 412.9, 412.9, -412.9],
                [348.3, 352.3, 352.3, 348.3, 348.3],
            ]
            SOLENOID_2 = [
                [-412.9, -412.9, 412.9, 412.9, -412.9],
                [364.9, 399.3, 399.3, 364.9, 364.9],
            ]
            SOLENOID_3 = [
                [-412.9, -412.9, 412.9, 412.9, -412.9],
                [425, 429, 429, 425, 425],
            ]

            MD1 = [
                [-563.8, -563.8, 563.8, 563.8, 417.9, 417.9, -417.9, -417.9, -563.8],
                [78.2, 645, 645, 78.2, 57.5, 446.1, 446.1, 57.5, 78.2],
            ]

            CONE1 = [
                [6.5, 230.7, 250.9, 412.9, 417.9, 563.8, 563.8, 6.5],
                [2.2, 31, 33.9, 56.8, 57.5, 78.2, 2.2, 2.2],
            ]
            CONE2 = [[-1 * i for i in CONE1[0]], CONE1[1]]

            BL = [[-563.8, -563.8, 563.8, 563.8, -563.8], [-2.2, 2.2, 2.2, -2.2, -2.2]]

            components = [
                ECAL1,
                ECAL2,
                ECAL3,
                HCAL1,
                HCAL2,
                HCAL3,
                SOLENOID,
                SOLENOID_2,
                SOLENOID_3,
                MD1,
                CONE1,
                CONE2,
                BL,
            ]
            cols = (
                ["dimgrey"] * 3
                + ["darkgrey"] * 3
                + ["gray"]
                + ["lightgrey"]
                + ["gray"]
                + ["grey"]
                + ["black"] * 2
                + ["white"]
            )

            for i, component in enumerate(components):
                # ax.plot(component[0], component[1], color=cols[i], lw=0.5)
                ax.fill_between(
                    component[0], component[1], color=cols[i], alpha=1, lw=0.1
                )
                new_y = [-1 * k for k in component[1]]
                # ax.plot(component[0], new_y, color=cols[i], lw=0.5)
                ax.fill_between(component[0], new_y, color=cols[i], alpha=1, lw=0.1)

            if xl:
                ax.set_xlabel("z-coordinate (cm)")

            if yl:
                ax.set_ylabel("y-coordinate (cm)")

            # ax.set_xlim(-564, 564)
            # ax.set_ylim(-645, 645)

    else:
        print("this geometry has not been implemented yet!")

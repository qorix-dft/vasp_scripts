#!/usr/bin/env python3
"""
Interactive per-layer 3D bar plot of the per-atom |dr|, matplotlib version.

Two-file variant working on CONTCAR geometries: an ab initio reference and an
MLFF-produced structure (also written from an ab initio-style output). For each
atom we take the minimum-image displacement between the two structures and plot
its magnitude

    dr = || r(CONTCAR_ai) - r(CONTCAR_mlff) ||   per atom

Bars sit at each atom's (x, y) Cartesian position with height = dr, one 3D
panel per atomic layer (trilayer slab). Opens a native draggable window
(left-drag rotate, scroll zoom).

The slab contains one interstitial atom (index hardcoded below): it is pulled
out of the layer panels and shown on its own so it does not distort the layer
splitting or the shared colour/height scales.

Bars are drawn as cylinders (circular footprint) instead of square prisms.

Run on WSL2 / Windows 11 (WSLg):
    pip install PyQt5
    QT_QPA_PLATFORM=wayland python3 delta_pos.py

Note: keep this script outside any directory that is rsync'd with --delete
from the supercomputer, or it will be removed on the next sync.
"""

import numpy as np
import matplotlib

matplotlib.use("QtAgg")

import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================
# User settings
# ==========================

CONTCAR_AI = "CONTCAR_ai"        # ab initio reference structure
CONTCAR_MLFF = "CONTCAR_mlff"    # MLFF-relaxed structure

# 1-based index of the interstitial atom as it appears in the CONTCAR ordering
# (i.e. the same numbering VASP / VESTA show). Set this to your interstitial.
INTERSTITIAL_INDEX = 1

# Geometry used to place the bars and to assign layers.
GEOMETRY_REF = CONTCAR_AI

N_LAYERS = 3

# Footprint radius of each cylindrical bar in the xy plane (Angstrom).
BAR_RADIUS = 0.22
BAR_NSIDES = 24

CMAP = "viridis"

# Bars with dr below this threshold (Ang) are drawn in LOW_COLOR instead of on
# the colormap, to flag atoms that barely move between the two structures.
LOW_DR_THRESHOLD = 0.02
LOW_COLOR = "red"

# Per-atom displacements are written here, one row per atom.
OUT_DAT = "dr_per_atom.dat"

# One window with all layers side by side (True), or one window per layer that
# you close to advance (False). The interstitial always gets its own window.
COMBINED_WINDOW = True

# ==========================


def read_contcar(filename):
    """Return Cartesian positions (natoms, 3) and the 3x3 lattice matrix (rows
    are lattice vectors) from a CONTCAR / POSCAR file."""
    with open(filename, "r", errors="ignore") as f:
        lines = f.readlines()

    scale = float(lines[1].split()[0])
    lattice = np.array([[float(x) for x in lines[i].split()[:3]]
                        for i in (2, 3, 4)], dtype=float)
    lattice *= scale

    # Line 5 is either element symbols (VASP5+) or the counts (VASP4).
    tok = lines[5].split()
    if all(t.lstrip("-").isdigit() for t in tok):
        counts_line = 5
    else:
        counts_line = 6
    counts = [int(x) for x in lines[counts_line].split()]
    natoms = sum(counts)

    idx = counts_line + 1
    if lines[idx].strip()[:1] in ("s", "S"):   # Selective dynamics
        idx += 1

    mode = lines[idx].strip()[:1].lower()
    cartesian = mode in ("c", "k")
    idx += 1

    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]]
                       for i in range(natoms)], dtype=float)

    if cartesian:
        cart = coords * scale
        frac = cart @ np.linalg.inv(lattice)
    else:
        frac = coords
        cart = frac @ lattice

    return cart, frac, lattice


def min_image_displacement(frac_a, frac_b, lattice):
    """Per-atom minimum-image displacement magnitude || r_a - r_b || (Ang),
    wrapping the fractional difference into [-0.5, 0.5) before going to
    Cartesian so atoms that cross a cell boundary are handled correctly."""
    dfrac = frac_a - frac_b
    dfrac -= np.round(dfrac)
    dcart = dfrac @ lattice
    return np.linalg.norm(dcart, axis=1), dcart


def assign_layers(z, n_layers):
    """Split atoms into n_layers groups by cutting the sorted z at the largest
    gaps. Returns integer labels (bottom to top) and each layer's mean z."""
    order = np.argsort(z)
    z_sorted = z[order]

    if n_layers < 2:
        return np.zeros(len(z), dtype=int), np.array([z.mean()])

    gaps = np.diff(z_sorted)
    cut_after = np.sort(np.argsort(gaps)[-(n_layers - 1):])

    labels_sorted = np.zeros(len(z), dtype=int)
    layer, prev = 0, 0
    for c in cut_after:
        labels_sorted[prev:c + 1] = layer
        layer += 1
        prev = c + 1
    labels_sorted[prev:] = layer

    labels = np.empty(len(z), dtype=int)
    labels[order] = labels_sorted
    z_means = np.array([z[labels == k].mean() for k in range(n_layers)])
    return labels, z_means


def write_per_atom(filename, positions, dr, dcart, labels, z_means,
                   interstitial0):
    """One row per atom: index, Cartesian position, layer, dx/dy/dz and dr."""
    natoms = len(positions)
    order = sorted(range(natoms),
                   key=lambda i: (labels[i] if i != interstitial0 else N_LAYERS, i))
    with open(filename, "w") as f:
        f.write("# Per-atom displacement dr = || r(CONTCAR_ai) - r(CONTCAR_mlff) ||\n")
        f.write("# (minimum-image convention; dx/dy/dz are Cartesian components)\n")
        f.write(f"# Interstitial atom index (1-based): {interstitial0 + 1}\n")
        f.write("# Layers by mean z (Ang): "
                + ", ".join(f"L{k+1}={z_means[k]:.4f}" for k in range(N_LAYERS))
                + "\n")
        f.write(f"# low_dr = 1 where dr < {LOW_DR_THRESHOLD:.3f} Ang "
                f"(drawn {LOW_COLOR})\n")
        f.write("#\n")
        f.write(f"#{'atom':>6} {'x':>12} {'y':>12} {'z':>12} "
                f"{'layer':>7} {'dx':>13} {'dy':>13} {'dz':>13} "
                f"{'dr':>16} {'low_dr':>7}\n")
        for i in order:
            layer_str = "inter" if i == interstitial0 else f"{labels[i] + 1}"
            f.write(f" {i + 1:>6d} "
                    f"{positions[i, 0]:12.6f} {positions[i, 1]:12.6f} "
                    f"{positions[i, 2]:12.6f} {layer_str:>7} "
                    f"{dcart[i, 0]:13.6f} {dcart[i, 1]:13.6f} {dcart[i, 2]:13.6f} "
                    f"{dr[i]:16.8e} {int(dr[i] < LOW_DR_THRESHOLD):7d}\n")


def draw_cylinders(ax, x, y, df, facecolors, radius, nsides):
    """Draw one vertical cylinder per atom: circular footprint at (x, y),
    height df. Faces = curved side (quads) + top cap (n-gon)."""
    theta = np.linspace(0.0, 2.0 * np.pi, nsides + 1)
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    side_polys, side_cols = [], []
    cap_polys, cap_cols = [], []
    for xi, yi, hi, col in zip(x, y, df, facecolors):
        cx = xi + radius * cos_t
        cy = yi + radius * sin_t
        for j in range(nsides):
            side_polys.append([(cx[j], cy[j], 0.0),
                               (cx[j + 1], cy[j + 1], 0.0),
                               (cx[j + 1], cy[j + 1], hi),
                               (cx[j], cy[j], hi)])
            side_cols.append(col)
        cap_polys.append([(cx[j], cy[j], hi) for j in range(nsides)])
        cap_cols.append(col)

    if side_polys:
        ax.add_collection3d(Poly3DCollection(
            side_polys, facecolors=side_cols, edgecolors=side_cols,
            linewidths=0.0, shade=True))
    if cap_polys:
        ax.add_collection3d(Poly3DCollection(
            cap_polys, facecolors=cap_cols,
            edgecolors="k", linewidths=0.2, shade=True))


def draw_layer(ax, x, y, df, norm, cmap, xlim, ylim, zmax, title):
    """One 3D cylinder chart: bars at (x, y), heights df."""
    facecolors = cmap(norm(df))
    facecolors[df < LOW_DR_THRESHOLD] = colors.to_rgba(LOW_COLOR)

    draw_cylinders(ax, x, y, df, facecolors, BAR_RADIUS, BAR_NSIDES)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(0.0, zmax)

    ax.set_xlabel(r"$x$  ($\mathrm{\AA}$)", labelpad=6)
    ax.set_ylabel(r"$y$  ($\mathrm{\AA}$)", labelpad=6)
    ax.set_zlabel(r"$|\Delta r|$  ($\mathrm{\AA}$)", labelpad=6)
    ax.set_title(title, fontsize=11)
    ax.view_init(elev=28, azim=-60)


def draw_interstitial(dr_i, dcart_i, norm, cmap, zmax):
    """Own window for the single interstitial atom: a circular bar for the
    total |dr| next to a small breakdown of the (dx, dy, dz) components."""
    fig = plt.figure(figsize=(9.6, 5.4))

    # Left: circular bar in the same style/scale as the layer panels.
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    col = (colors.to_rgba(LOW_COLOR) if dr_i < LOW_DR_THRESHOLD
           else cmap(norm(dr_i)))
    draw_cylinders(ax, np.array([0.0]), np.array([0.0]), np.array([dr_i]),
                   [col], BAR_RADIUS * 3.0, BAR_NSIDES)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_zlim(0.0, zmax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zlabel(r"$|\Delta r|$  ($\mathrm{\AA}$)", labelpad=6)
    ax.set_title(rf"interstitial   $|\Delta r|$ = {dr_i:.4f} $\mathrm{{\AA}}$",
                 fontsize=11)
    ax.view_init(elev=18, azim=-60)

    # Right: signed Cartesian components, so you can see direction of the shift.
    ax2 = fig.add_subplot(1, 2, 2)
    comp = ["$\\Delta x$", "$\\Delta y$", "$\\Delta z$"]
    vals = dcart_i
    bar_cols = ["#4c72b0", "#dd8452", "#55a868"]
    ax2.bar(comp, vals, color=bar_cols, edgecolor="k", linewidth=0.4, width=0.6)
    ax2.axhline(0.0, color="k", lw=0.8)
    lim = max(0.02, np.abs(vals).max() * 1.3)
    ax2.set_ylim(-lim, lim)
    ax2.set_ylabel(r"displacement  ($\mathrm{\AA}$)")
    ax2.set_title("component breakdown", fontsize=11)
    for xi, v in enumerate(vals):
        ax2.text(xi, v + np.sign(v) * lim * 0.03, f"{v:+.3f}",
                 ha="center", va="bottom" if v >= 0 else "top", fontsize=9)

    fig.suptitle("Interstitial displacement", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))


def main():
    cart_ai, frac_ai, lat_ai = read_contcar(CONTCAR_AI)
    cart_ml, frac_ml, lat_ml = read_contcar(CONTCAR_MLFF)

    natoms = len(cart_ai)
    if len(cart_ml) != natoms:
        raise ValueError(
            f"Atom count mismatch: {CONTCAR_AI} has {natoms}, "
            f"{CONTCAR_MLFF} has {len(cart_ml)}"
        )

    interstitial0 = INTERSTITIAL_INDEX - 1
    if not (0 <= interstitial0 < natoms):
        raise ValueError(
            f"INTERSTITIAL_INDEX={INTERSTITIAL_INDEX} out of range "
            f"(1..{natoms})"
        )

    positions = cart_ai if GEOMETRY_REF == CONTCAR_AI else cart_ml
    lattice = lat_ai if GEOMETRY_REF == CONTCAR_AI else lat_ml
    frac_ref = frac_ai if GEOMETRY_REF == CONTCAR_AI else frac_ml

    dr, dcart = min_image_displacement(frac_ai, frac_ml, lattice)

    # Layer assignment on the framework atoms only (drop the interstitial so it
    # doesn't create a spurious z-gap).
    frame = np.ones(natoms, dtype=bool)
    frame[interstitial0] = False
    labels_frame, z_means = assign_layers(positions[frame, 2], N_LAYERS)
    labels = np.full(natoms, -1, dtype=int)
    labels[frame] = labels_frame

    # Shared height/colour scale across all layers AND the interstitial.
    dmax = float(dr.max()) if dr.max() > 0 else 1.0
    zmax = dmax * 1.05
    norm = colors.Normalize(vmin=0.0, vmax=dmax)
    cmap = plt.get_cmap(CMAP)

    pad = 1.0
    fx = positions[frame]
    xlim = (fx[:, 0].min() - pad, fx[:, 0].max() + pad)
    ylim = (fx[:, 1].min() - pad, fx[:, 1].max() + pad)

    write_per_atom(OUT_DAT, positions, dr, dcart, labels, z_means, interstitial0)

    n_low = int((dr[frame] < LOW_DR_THRESHOLD).sum())
    print(f"Atoms: {natoms} (framework {int(frame.sum())} + 1 interstitial) | "
          f"max |dr| = {dr.max():.6f} A (atom {int(dr.argmax())+1})")
    print(f"Interstitial (atom {INTERSTITIAL_INDEX}): |dr| = {dr[interstitial0]:.6f} A, "
          f"(dx,dy,dz)=({dcart[interstitial0,0]:+.4f}, "
          f"{dcart[interstitial0,1]:+.4f}, {dcart[interstitial0,2]:+.4f}) A")
    for k in range(N_LAYERS):
        mk = labels == k
        print(f"  layer {k+1}: {int(mk.sum()):4d} atoms, <z>={z_means[k]:.4f} A, "
              f"sum dr={dr[mk].sum():.6f}, max dr={dr[mk].max():.6f}")
    print(f"  framework atoms with dr < {LOW_DR_THRESHOLD:.3f} A: {n_low} "
          f"(drawn {LOW_COLOR})")
    print(f"Wrote {OUT_DAT}")

    def title_for(k):
        m = labels == k
        return (f"layer {k + 1}   "
                rf"$\langle z \rangle$ = {z_means[k]:.2f} $\mathrm{{\AA}}$"
                f"   ({int(m.sum())} atoms)")

    if COMBINED_WINDOW:
        fig = plt.figure(figsize=(6.2 * N_LAYERS, 6.4))
        for k in range(N_LAYERS):
            m = labels == k
            ax = fig.add_subplot(1, N_LAYERS, k + 1, projection="3d")
            draw_layer(ax, positions[m, 0], positions[m, 1], dr[m],
                       norm, cmap, xlim, ylim, zmax, title_for(k))

        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(mappable, ax=fig.axes, shrink=0.65, pad=0.04)
        cbar.set_label(r"$|\Delta r|$  ($\mathrm{\AA}$)")

        fig.suptitle(r"Per-atom $|\Delta r|$ by layer "
                     r"(cylinders at atomic $(x, y)$ positions)", y=0.98)
    else:
        for k in range(N_LAYERS):
            m = labels == k
            fig = plt.figure(figsize=(7.4, 6.6))
            ax = fig.add_subplot(111, projection="3d")
            draw_layer(ax, positions[m, 0], positions[m, 1], dr[m],
                       norm, cmap, xlim, ylim, zmax, title_for(k))

            mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.10)
            cbar.set_label(r"$|\Delta r|$  ($\mathrm{\AA}$)")

    draw_interstitial(dr[interstitial0], dcart[interstitial0], norm, cmap, zmax)
    plt.show()


if __name__ == "__main__":
    main()

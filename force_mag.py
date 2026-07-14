#!/usr/bin/env python3
"""
Interactive per-layer 3D bar plot of the per-atom |F|, matplotlib version.

Single-file variant: one OUTCAR, and the per-atom force magnitude

    F = || F(OUTCAR) ||   per atom

taken from the final POSITION / TOTAL-FORCE block. Bars sit at each atom's
(x, y) Cartesian position with height = |F|, one 3D panel per atomic layer
(trilayer slab). Opens a native draggable window (left-drag rotate, scroll
zoom).

The slab contains one interstitial atom (index hardcoded below): it is pulled
out of the layer panels and shown on its own so it does not distort the layer
splitting or the shared colour/height scales.

Bars are drawn as cylinders (circular footprint) instead of square prisms.

Run on WSL2 / Windows 11 (WSLg):
    pip install PyQt5
    QT_QPA_PLATFORM=wayland python3 force_mag.py

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

OUTCAR = "OUTCAR"

# 1-based index of the interstitial atom as it appears in the OUTCAR ordering
# (i.e. the same numbering VASP / VESTA show). Set this to your interstitial.
INTERSTITIAL_INDEX = 1

N_LAYERS = 3

# Footprint radius of each cylindrical bar in the xy plane (Angstrom).
BAR_RADIUS = 0.22
BAR_NSIDES = 24

CMAP = "viridis"

# Bars with |F| below this threshold (eV/Ang) are drawn in LOW_COLOR instead of
# on the colormap, to flag atoms that are essentially force-free (converged).
LOW_F_THRESHOLD = 0.05
LOW_COLOR = "red"

# Per-atom forces are written here, one row per atom.
OUT_DAT = "F_per_atom.dat"

# One window with all layers side by side (True), or one window per layer that
# you close to advance (False). The interstitial always gets its own window.
COMBINED_WINDOW = True

# ==========================


def read_last_force_block(filename):
    """Return (positions, forces), each (natoms, 3), from the final
    POSITION / TOTAL-FORCE block of an OUTCAR."""
    with open(filename, "r", errors="ignore") as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if "POSITION" in line and "TOTAL-FORCE" in line:
            header_idx = i

    if header_idx is None:
        raise ValueError(f"No POSITION/TOTAL-FORCE block found in {filename}")

    rows = []
    started = False
    for line in lines[header_idx + 1:]:
        parts = line.split()
        try:
            values = [float(x) for x in parts]
        except ValueError:
            values = None

        if values is not None and len(values) == 6:
            rows.append(values)
            started = True
        elif started:
            break

    if not rows:
        raise ValueError(f"Force block in {filename} contained no data rows")

    block = np.array(rows, dtype=float)
    return block[:, 0:3], block[:, 3:6]


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


def write_per_atom(filename, positions, F, forces, labels, z_means,
                   interstitial0):
    """One row per atom: index, Cartesian position, layer, Fx/Fy/Fz and |F|."""
    natoms = len(positions)
    order = sorted(range(natoms),
                   key=lambda i: (labels[i] if i != interstitial0 else N_LAYERS, i))
    with open(filename, "w") as f:
        f.write("# Per-atom force magnitude F = || F(OUTCAR) ||\n")
        f.write("# (Fx/Fy/Fz are the Cartesian force components, eV/Ang)\n")
        f.write(f"# Interstitial atom index (1-based): {interstitial0 + 1}\n")
        f.write("# Layers by mean z (Ang): "
                + ", ".join(f"L{k+1}={z_means[k]:.4f}" for k in range(N_LAYERS))
                + "\n")
        f.write(f"# low_F = 1 where |F| < {LOW_F_THRESHOLD:.3f} eV/Ang "
                f"(drawn {LOW_COLOR})\n")
        f.write("#\n")
        f.write(f"#{'atom':>6} {'x':>12} {'y':>12} {'z':>12} "
                f"{'layer':>7} {'Fx':>13} {'Fy':>13} {'Fz':>13} "
                f"{'|F|':>16} {'low_F':>6}\n")
        for i in order:
            layer_str = "inter" if i == interstitial0 else f"{labels[i] + 1}"
            f.write(f" {i + 1:>6d} "
                    f"{positions[i, 0]:12.6f} {positions[i, 1]:12.6f} "
                    f"{positions[i, 2]:12.6f} {layer_str:>7} "
                    f"{forces[i, 0]:13.6f} {forces[i, 1]:13.6f} {forces[i, 2]:13.6f} "
                    f"{F[i]:16.8e} {int(F[i] < LOW_F_THRESHOLD):6d}\n")


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
    facecolors[df < LOW_F_THRESHOLD] = colors.to_rgba(LOW_COLOR)

    draw_cylinders(ax, x, y, df, facecolors, BAR_RADIUS, BAR_NSIDES)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(0.0, zmax)

    ax.set_xlabel(r"$x$  ($\mathrm{\AA}$)", labelpad=6)
    ax.set_ylabel(r"$y$  ($\mathrm{\AA}$)", labelpad=6)
    ax.set_zlabel(r"$|F|$  (eV/$\mathrm{\AA}$)", labelpad=6)
    ax.set_title(title, fontsize=11)
    ax.view_init(elev=28, azim=-60)


def draw_interstitial(F_i, force_i, norm, cmap, zmax):
    """Own window for the single interstitial atom: a circular bar for the
    total |F| next to a small breakdown of the (Fx, Fy, Fz) components."""
    fig = plt.figure(figsize=(9.6, 5.4))

    # Left: circular bar in the same style/scale as the layer panels.
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    col = (colors.to_rgba(LOW_COLOR) if F_i < LOW_F_THRESHOLD
           else cmap(norm(F_i)))
    draw_cylinders(ax, np.array([0.0]), np.array([0.0]), np.array([F_i]),
                   [col], BAR_RADIUS * 3.0, BAR_NSIDES)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_zlim(0.0, zmax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zlabel(r"$|F|$  (eV/$\mathrm{\AA}$)", labelpad=6)
    ax.set_title(rf"interstitial   $|F|$ = {F_i:.4f} eV/$\mathrm{{\AA}}$",
                 fontsize=11)
    ax.view_init(elev=18, azim=-60)

    # Right: signed Cartesian components, so you can see the direction of the
    # residual force.
    ax2 = fig.add_subplot(1, 2, 2)
    comp = ["$F_x$", "$F_y$", "$F_z$"]
    vals = force_i
    bar_cols = ["#4c72b0", "#dd8452", "#55a868"]
    ax2.bar(comp, vals, color=bar_cols, edgecolor="k", linewidth=0.4, width=0.6)
    ax2.axhline(0.0, color="k", lw=0.8)
    lim = max(0.02, np.abs(vals).max() * 1.3)
    ax2.set_ylim(-lim, lim)
    ax2.set_ylabel(r"force  (eV/$\mathrm{\AA}$)")
    ax2.set_title("component breakdown", fontsize=11)
    for xi, v in enumerate(vals):
        ax2.text(xi, v + np.sign(v) * lim * 0.03, f"{v:+.3f}",
                 ha="center", va="bottom" if v >= 0 else "top", fontsize=9)

    fig.suptitle("Interstitial force", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))


def main():
    positions, forces = read_last_force_block(OUTCAR)
    natoms = len(positions)

    interstitial0 = INTERSTITIAL_INDEX - 1
    if not (0 <= interstitial0 < natoms):
        raise ValueError(
            f"INTERSTITIAL_INDEX={INTERSTITIAL_INDEX} out of range "
            f"(1..{natoms})"
        )

    F = np.linalg.norm(forces, axis=1)

    # Layer assignment on the framework atoms only (drop the interstitial so it
    # doesn't create a spurious z-gap).
    frame = np.ones(natoms, dtype=bool)
    frame[interstitial0] = False
    labels_frame, z_means = assign_layers(positions[frame, 2], N_LAYERS)
    labels = np.full(natoms, -1, dtype=int)
    labels[frame] = labels_frame

    # Shared height/colour scale across all layers AND the interstitial.
    fmax = float(F.max()) if F.max() > 0 else 1.0
    zmax = fmax * 1.05
    norm = colors.Normalize(vmin=0.0, vmax=fmax)
    cmap = plt.get_cmap(CMAP)

    pad = 1.0
    fx = positions[frame]
    xlim = (fx[:, 0].min() - pad, fx[:, 0].max() + pad)
    ylim = (fx[:, 1].min() - pad, fx[:, 1].max() + pad)

    write_per_atom(OUT_DAT, positions, F, forces, labels, z_means, interstitial0)

    n_low = int((F[frame] < LOW_F_THRESHOLD).sum())
    print(f"Atoms: {natoms} (framework {int(frame.sum())} + 1 interstitial) | "
          f"max |F| = {F.max():.6f} eV/A (atom {int(F.argmax())+1})")
    print(f"Interstitial (atom {INTERSTITIAL_INDEX}): |F| = {F[interstitial0]:.6f} eV/A, "
          f"(Fx,Fy,Fz)=({forces[interstitial0,0]:+.4f}, "
          f"{forces[interstitial0,1]:+.4f}, {forces[interstitial0,2]:+.4f}) eV/A")
    for k in range(N_LAYERS):
        mk = labels == k
        print(f"  layer {k+1}: {int(mk.sum()):4d} atoms, <z>={z_means[k]:.4f} A, "
              f"sum |F|={F[mk].sum():.6f}, max |F|={F[mk].max():.6f}")
    print(f"  framework atoms with |F| < {LOW_F_THRESHOLD:.3f} eV/A: {n_low} "
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
            draw_layer(ax, positions[m, 0], positions[m, 1], F[m],
                       norm, cmap, xlim, ylim, zmax, title_for(k))

        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(mappable, ax=fig.axes, shrink=0.65, pad=0.04)
        cbar.set_label(r"$|F|$  (eV/$\mathrm{\AA}$)")

        fig.suptitle(r"Per-atom $|F|$ by layer "
                     r"(cylinders at atomic $(x, y)$ positions)", y=0.98)
    else:
        for k in range(N_LAYERS):
            m = labels == k
            fig = plt.figure(figsize=(7.4, 6.6))
            ax = fig.add_subplot(111, projection="3d")
            draw_layer(ax, positions[m, 0], positions[m, 1], F[m],
                       norm, cmap, xlim, ylim, zmax, title_for(k))

            mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.10)
            cbar.set_label(r"$|F|$  (eV/$\mathrm{\AA}$)")

    draw_interstitial(F[interstitial0], forces[interstitial0], norm, cmap, zmax)
    plt.show()


if __name__ == "__main__":
    main()

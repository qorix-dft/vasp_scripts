ML Accuracy: Hexagonal boron nitride

# Per-Atom Displacement Viewer (`delta_pos.py`)

Interactive 3D visualisation of the per-atom displacement between two relaxed
structures — an **ab initio** reference and an **MLFF**-relaxed geometry — for a
trilayer slab containing a single interstitial atom.

For every atom the script computes the minimum-image displacement

```
dr = || r(CONTCAR_ai) - r(CONTCAR_mlff) ||
```

and draws it as a vertical cylinder standing at that atom's `(x, y)` position,
with height equal to `dr`. One 3D panel is produced per atomic layer, plus a
separate window for the interstitial.

## What it's for

Comparing an MLFF-relaxed structure against the ab initio reference tells you
*where* the force field is accurate and where it isn't. Plotting `dr` spatially
(rather than as a single RMSD number) makes it obvious whether the error is
uniformly small, localised around the defect, or concentrated in a particular
layer.

## Requirements

- Python 3
- `numpy`, `matplotlib`
- A Qt binding for the interactive window (`pip install PyQt5`)

The script sets `matplotlib.use("QtAgg")`, so it opens a native, draggable
window rather than writing a static image.

## Input

Two VASP structure files in the working directory:

| File           | Description                                    |
|----------------|------------------------------------------------|
| `CONTCAR_ai`   | Ab initio reference structure                  |
| `CONTCAR_mlff` | MLFF-relaxed structure                         |

Both must contain the **same number of atoms in the same order**; the script
raises an error on a count mismatch. Standard `CONTCAR`/`POSCAR` format is
handled, including:

- VASP 4 and VASP 5+ headers (element-symbol line present or absent)
- `Selective dynamics` lines
- Both `Direct`/fractional and `Cartesian` coordinates
- A global scaling factor on line 2

## Usage

```bash
QT_QPA_PLATFORM=wayland python3 delta_pos.py
```

The `QT_QPA_PLATFORM=wayland` prefix is for WSL2 / Windows 11 with WSLg. On a
native Linux desktop you can usually just run `python3 delta_pos.py`.

Interaction: **left-drag** to rotate, **scroll** to zoom.

> **Note:** keep this script outside any directory that is `rsync`'d from the
> supercomputer with `--delete`, or it will be removed on the next sync.

## Output

**Windows**

- One 3D panel per layer (all layers side by side in a single window by
  default), sharing a common height and colour scale.
- A separate window for the interstitial atom, showing its total `|dr|` as a
  cylinder alongside a signed `(dx, dy, dz)` component breakdown — so you can
  see the *direction* of its shift, not just the magnitude.

**File — `dr_per_atom.dat`**

One row per atom, sorted by layer (interstitial last), with columns:

```
atom   x   y   z   layer   dx   dy   dz   dr   low_dr
```

Header comments record the interstitial index, the mean z of each layer, and
the low-`dr` threshold used.

## Configuration

All settings are constants near the top of the script.

| Variable              | Default          | Description                                                     |
|-----------------------|------------------|-----------------------------------------------------------------|
| `CONTCAR_AI`          | `"CONTCAR_ai"`   | Ab initio reference filename                                     |
| `CONTCAR_MLFF`        | `"CONTCAR_mlff"` | MLFF structure filename                                          |
| `INTERSTITIAL_INDEX`  | `1`              | **1-based** index of the interstitial, in CONTCAR/VESTA numbering |
| `GEOMETRY_REF`        | `CONTCAR_AI`     | Which structure supplies the bar positions and layer assignment   |
| `N_LAYERS`            | `3`              | Number of layers to split the slab into                          |
| `BAR_RADIUS`          | `0.22`           | Cylinder footprint radius (Å)                                    |
| `BAR_NSIDES`          | `24`             | Polygon sides approximating each cylinder                        |
| `CMAP`                | `"viridis"`      | Colormap for `dr`                                                |
| `LOW_DR_THRESHOLD`    | `0.02`           | Atoms below this `dr` (Å) are flagged                            |
| `LOW_COLOR`           | `"red"`          | Colour used for those flagged atoms                              |
| `OUT_DAT`             | `"dr_per_atom.dat"` | Per-atom output filename                                      |
| `COMBINED_WINDOW`     | `True`           | All layers in one window, vs. one window per layer               |

### Notes on specific settings

- **`INTERSTITIAL_INDEX` is 1-based** — use the same number VASP or VESTA
  shows. It's converted internally.
- **`GEOMETRY_REF`** is compared by string value against `CONTCAR_AI`, so set it
  to one of the two filename constants (not an arbitrary path).
- **`LOW_DR_THRESHOLD`** only changes how bars are *coloured* — those atoms are
  still drawn at their true height and still contribute to the shared scale.

## How it works

1. **Parsing** — both structures are read into Cartesian and fractional
   coordinates plus the lattice matrix.
2. **Displacement** — the fractional difference is wrapped into `[-0.5, 0.5)`
   before converting to Cartesian, so atoms that drift across a periodic
   boundary give the correct small displacement rather than a full cell vector.
3. **Layer assignment** — atoms are sorted by `z` and cut at the `N_LAYERS - 1`
   largest gaps. The interstitial is **excluded** from this step, since sitting
   between layers it would otherwise open a spurious gap and corrupt the
   splitting.
4. **Shared scaling** — height limits and the colour normalisation are computed
   once across *all* atoms including the interstitial, so every panel and the
   interstitial window are directly comparable.
5. **Rendering** — each bar is built as an explicit `Poly3DCollection`: a curved
   side surface of quads plus a top cap. The interstitial is drawn at 3× radius
   so it's visible on its own.

## Console summary

On each run the script prints atom counts, the global maximum `|dr|` and which
atom it belongs to, the interstitial's displacement and components, and then
per layer: atom count, mean `z`, summed `dr` and maximum `dr`.

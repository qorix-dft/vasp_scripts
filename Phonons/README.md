# Phonons

Lightweight sanity-checking tools for [Phonopy](https://phonopy.github.io/phonopy/) output.

## Overview

`compare_band_yaml.py` compares two Phonopy `band.yaml` files — typically a small
reference cell against a larger supercell of the same material — and reports
whether they describe consistent structures and consistent phonon spectra.

The idea is to catch common failures **before** committing to expensive
eigenvector-overlap or unfolding analysis. It deliberately parses only:

- metadata (`natom`, lattice vectors, q-points)
- atomic species and fractional coordinates
- mode frequencies

No eigenvectors are read, so the script stays fast even on large `band.yaml`
files.

## What it checks

### 1. Metadata summary
For each file: atom count, expected number of modes per q-point (`3N`), number
of q-points, species composition, lattice vector lengths, and the q-point list.

### 2. Frequency summary
Per q-point, prints min / median / max frequency plus counts of:

- `negative < -0.01 THz` — soft or imaginary modes
- `negative < -0.1 THz` — strongly imaginary modes (real instabilities vs. noise)
- `|freq| <= 0.1 THz` — acoustic / near-zero modes

The ten lowest frequencies are printed so acoustic sum-rule violations are
immediately visible.

### 3. Position embedding
Checks that every atom in the small cell has a species-matched counterpart in
the large cell. Matching is greedy, one-to-one, and uses the **minimum-image
convention** in Cartesian space (via the large cell's lattice), so periodic
wrap-around does not produce false mismatches.

Reports matched/unmatched counts, unused large-cell atoms, and the
min/mean/p95/max matching distance. A large `max` distance or any unmatched
atoms means the two cells are not the structures you think they are.

### 4. Frequency nearest-neighbour check
For each positive mode in the small cell, finds the closest positive mode in the
large cell and reports the distribution of those differences
(mean, p50, p90, p99, max). If the small cell's spectrum is genuinely a subset
of the large cell's, these differences should be near zero.

## Requirements

- Python ≥ 3.9
- NumPy

```bash
pip install numpy
```

## Usage

```bash
python compare_band_yaml.py SMALL_BAND_YAML LARGE_BAND_YAML [--position-tolerance TOL]
```

**Arguments**

| Argument | Description |
|---|---|
| `small_band_yaml` | Path to the reference / smaller cell `band.yaml` |
| `large_band_yaml` | Path to the supercell / larger cell `band.yaml` |
| `--position-tolerance` | Max distance (Å) for an atom match. Default `0.15` |

**Example**

```bash
python compare_band_yaml.py unitcell/band.yaml supercell_2x2x1/band.yaml
python compare_band_yaml.py unitcell/band.yaml supercell_2x2x1/band.yaml --position-tolerance 0.05
```

## Interpreting the output

| Symptom | Likely cause |
|---|---|
| Unmatched atoms, large max distance | Wrong supercell, different relaxation, or origin shift |
| Many unused large-cell atoms | Expected — the supercell has more atoms than the reference |
| Large `negative < -0.1` count | Structural instability, or an unconverged force constant calculation |
| Near-zero modes ≠ 3 at Γ | Acoustic sum rule not applied / broken |
| Large p99 or max frequency difference | Spectra do not correspond; check q-point mapping and cell commensurability |

## Notes

- Frequencies are assumed to be in THz, as written by Phonopy.
- The parser is regex-based rather than a full YAML load; this is intentional
  for speed on multi-MB `band.yaml` files, but it assumes standard Phonopy
  formatting.
- Atom matching is greedy and order-dependent. For structures with atoms much
  closer together than `--position-tolerance`, tighten the tolerance.


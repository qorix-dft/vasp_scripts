 Phonons

Lightweight sanity-checking tools for [Phonopy](https://phonopy.github.io/phonopy/) output.


# 'compare_band_yaml.py'
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

## Theory: the eigenvector correlation function

The intended follow-up analysis (once the cheap checks above pass) is a pairwise
comparison of phonon eigenvectors between two structures — for example, a
high-symmetry planar reference and a symmetry-broken, out-of-plane distorted
structure. Modes are compared using the correlation function

$$
\mathrm{cor}(k, k') = \left| \sum_{i,\alpha} \Delta r_{k;\alpha i}\, \Delta r_{k';\alpha i} \right|
$$

where $\Delta r_{k;\alpha i}$ is the component of the normalized eigenvector of
the $k$-th vibrational mode for atom $\alpha$ along direction
$i \in \{x, y, z\}$.

Correlation values run from 0 to 1: **0** means the modes $k$ and $k'$ are
orthogonal, **1** means they are identical. Interpretation:

- Non-zero elements concentrated on the **diagonal** with values near 1 → the
  mode is essentially unchanged between the two structures (expected for bulk
  phonons, which are only weakly perturbed).
- **Off-diagonal** weight → the mode's polarization changes between structures.
  A single mode in one structure mapping onto several modes in the other means
  it has *split* rather than merely shifted in energy — the characteristic
  signature of a soft mode decomposing into multiple quasi-local modes upon
  symmetry-breaking distortion.

A related useful quantity is the participation ratio, which quantifies how
localized a mode is:

$$
\mathrm{PR}_k = \sum_{\alpha} \left( \sum_{i} \Delta r_{k;\alpha i}^{2} \right)^{2}
$$

High PR indicates that only a few atoms participate appreciably — i.e. a
quasi-local mode rather than a delocalized bulk phonon.

> **Note:** the current script does not read eigenvectors, so it does not yet
> compute `cor(k, k')` or `PR_k`. It validates that the two `band.yaml` files
> are comparable in the first place, which is the prerequisite for this analysis.

### Reference

The correlation function and participation ratio above follow the formulation
in:

> Z. Tang, F. Jia, G. J. Kim-Reyes, Y. Wu, J. R. Chelikowsky, and P. Zhang,
> *Out-of-plane displacement of quantum color centers in monolayer h-BN*,
> arXiv:2503.06931 [cond-mat.mes-hall] (submitted 10 March 2025).
> https://arxiv.org/abs/2503.06931

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


# VASP Ab Initio Run Script Generator

`generate_run_scripts.sh` splits a batch of ab initio VASP relaxations (one
per `POSCAR-*` displacement) into groups, and writes a ready-to-submit
scheduler script for each group. It's scheduler-agnostic — SLURM, PBS/Torque,
LSF, or anything else — via a single editable settings block.

## Why

VASP jobs are often limited by a scheduler's max walltime. If you have
hundreds or thousands of displacement configurations to relax, one giant job
won't finish in time. This script divides the work into fixed-size groups,
each becoming its own submittable job, so you can process the whole set
across many separate allocations.

## Requirements

- Bash
- A directory containing:
  - `POSCAR-1`, `POSCAR-2`, ... `POSCAR-N` (zero-padded to 3 digits for
    indices under 100 — e.g. `POSCAR-007`, `POSCAR-042`, `POSCAR-1080`)
  - `KPOINTS`, `POTCAR`, `INCAR`
- A working scheduler command on your cluster (`sbatch`, `qsub`, `bsub`, ...)
  and VASP already set up (module/binary path)

## Usage

```bash
bash generate_run_scripts.sh              # generate, then asks y/n to submit
bash generate_run_scripts.sh --submit     # generate + submit everything, no prompt
bash generate_run_scripts.sh --no-submit  # generate only, never submit
```

Manual submission always works too:

```bash
sbatch run_01.sh   # or qsub / bsub / whatever SUBMIT_CMD is set to
```

## What it creates

For a run covering displacements 1–60, `run_01.sh` will contain a loop that,
when executed by the scheduler, for each displacement `i`:

1. Creates a `disp-<i>/` directory (e.g. `disp-001`, `disp-060`)
2. Moves `POSCAR-<i>` into it as `POSCAR`
3. Copies `KPOINTS`, `POTCAR`, `INCAR` into it
4. Launches VASP inside that directory, writing `vasp.out`

Each generated `run_XX.sh` is a complete, independent job — there's no
cross-job dependency, so they can all be queued and run concurrently.

> **Note:** `POSCAR-<i>` files are *moved*, not copied, into their
> displacement directories — the originals are consumed once a group runs.

## Configuration

All settings are edited directly in the script. There are three blocks:

### 1. Submission behavior

```bash
SUBMIT_DEFAULT="ask"   # "ask" | "yes" | "no"
```

Controls what happens when the script is run with no `--submit`/`--no-submit`
flag.

### 2. Supercomputer / scheduler settings

The block to replace when moving to a different machine:

```bash
SUBMIT_CMD="sbatch"          # sbatch, qsub, bsub, ...
JOB_NAME_PREFIX="phabinit"   # group number is appended automatically (g01, g02, ...)

JOB_HEADER=$(cat << 'HEADER_EOF'
#!/bin/bash -l
#SBATCH --job-name=__JOBNAME__
#SBATCH --account=...
...
HEADER_EOF
)

LAUNCH_CMD=$(cat << 'LAUNCH_EOF'
  srun ... "$VASP" > vasp.out
LAUNCH_EOF
)
```

- `JOB_HEADER` is copied verbatim into every generated run file — scheduler
  directives, module loads, and environment exports all go here. Use the
  placeholder `__JOBNAME__` wherever the per-group job name should be
  inserted.
- `LAUNCH_CMD` is the actual parallel-launch line that runs VASP
  (`srun`, `mpirun`, `aprun`, `ibrun`, ...).
- Both use a quoted heredoc (`<< 'HEADER_EOF'`), so nothing inside is
  evaluated when the generator runs — write runtime variables
  (`$SLURM_SUBMIT_DIR`, `$PBS_O_WORKDIR`, `$VASP`, etc.) exactly as they
  should appear; they're only evaluated later, when the generated job runs
  on the cluster.

The script is pre-filled with Pawsey/SLURM settings, so it works as-is on
Pawsey. Swap this block out entirely to target a different cluster.

### 3. Grouping parameters

```bash
TOTAL_POSCARS=1080   # total number of POSCAR-* files to process
GROUP_SIZE=60         # displacements handled per run file/job
OUTDIR="."            # where generated run_XX.sh files are written
```

`GROUP_SIZE` should be chosen so a group's total VASP runtime fits inside
your scheduler's walltime limit. `TOTAL_POSCARS` and `GROUP_SIZE` together
determine how many `run_XX.sh` files get created.

## How grouping works

Given `TOTAL_POSCARS` and `GROUP_SIZE`, the number of groups is:

```
NUM_GROUPS = ceil(TOTAL_POSCARS / GROUP_SIZE)
```

Groups are contiguous and non-overlapping (e.g. with `GROUP_SIZE=60`:
`1–60`, `61–120`, ..., with the final group truncated to `TOTAL_POSCARS` if
it doesn't divide evenly).

## Files

| File                     | Description                                      |
|--------------------------|---------------------------------------------------|
| `generate_run_scripts.sh`| The generator (this script)                       |
| `run_01.sh` ... `run_NN.sh` | Generated, per-group scheduler scripts (output) |

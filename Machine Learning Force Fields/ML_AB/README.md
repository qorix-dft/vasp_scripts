# Machine Learning Force Fields: ML_AB

Tools for VASP machine-learned force field (MLFF) training data.

# rm_structures.f90

Removes a section (or several sections) of configurations from an `ML_ABN`
training-data file and renumbers the remaining structures and basis sets
accordingly. Useful for throwing away a bad segment of an MD training run,
e.g. structures 30 to 100 out of 3000.

```
gfortran -o purge_structures purge_structures.f90
./purge_structures
```

The program reads `ML_ABN` from the current directory, prints the total
number of configurations, and prompts for the first and last configuration
of the section to eliminate (e.g. `30 100`). Sections can be entered
repeatedly; enter `0 0` when done. The purged data set is written to
`ML_ABN_purge` with:

- the eliminated configurations removed and the remaining ones renumbered
  consecutively (`Configuration num.`),
- basis sets pointing to eliminated structures removed, and the remaining
  entries remapped to the new configuration numbers,
- the number of configurations (line 5) and the numbers of basis sets per
  atom type updated automatically — no manual editing needed.

Copy/rename the `ML_ABN*` files as needed afterwards.



# make_ab_outcars

Build a VASP `ML_AB` training file (for `ML_MODE = select`) from a series of
`OUTCAR` files.

## What it does

- Scans `CONF/OUTCAR1`, `CONF/OUTCAR2`, ... (consecutive numbering, stops at
  the first missing file) and detects the atom types, number of atom types,
  atoms per type, and atomic masses **directly from each OUTCAR** — nothing
  is hardcoded. Different OUTCARs may contain different systems/compositions;
  the `ML_AB` header gets the union of all atom types found, and each
  configuration keeps the types/atom-counts of its own OUTCAR.
- Writes configurations to a scratch file first and copies them to `ML_AB`
  once the final count is known, so the configuration count on line 5 of
  `ML_AB` is always correct — no manual editing needed.
- From each OUTCAR, keeps only the **last structure** plus `NBACK` structures
  further back (`NBACK = 0` → last structure only, `NBACK = 2` → last 3 ionic
  steps, etc.). If an OUTCAR has fewer steps than that, all of them are kept.
- Assigns a system name to each configuration in one of two ways:
  1. **From a names file** (default `CONF/names.dat`): line *i* is the name
     for `CONF/OUTCARi` (blank lines and lines starting with `#` are
     skipped; names may contain blanks).
  2. **From ranges entered interactively**, e.g. `1 100 graphite`, one range
     per line, terminated with `0 0 end` (names entered this way must not
     contain blanks).
- Configurations are written to `ML_AB` **grouped by name**, in the order the
  names first appear.
- Also writes `energy.dat`, a plain list of the total energy (eV) of every
  configuration written to `ML_AB`, in the same order.

## Requirements

- `gfortran` (no external libraries)

## Compilation

```bash
gfortran -o make_ab_outcars make_ab_outcars.f90
```

## Input layout

```
CONF/
├── OUTCAR1
├── OUTCAR2
├── OUTCAR3
├── ...
└── names.dat        # optional, only needed for naming mode 1
```

OUTCAR numbering must start at 1 and be consecutive — the program stops
scanning at the first missing `CONF/OUTCARi`.

> If you collected your OUTCARs with a helper script that names them
> `outcars/OUTCAR1`, `outcars/OUTCAR2`, ..., just rename/symlink that folder
> to `CONF` before running this program.

`names.dat` format (naming mode 1):

```
# lines starting with '#' and blank lines are skipped
graphite
graphite
bulk_hbn
bulk_hbn
```

Line *i* gives the name for `CONF/OUTCARi`.

## Usage

```bash
./make_ab_outcars
```

The program will prompt for:

1. **NBACK** — how many extra structures to take back from the last one
   (`0` = last structure only).
2. **Naming mode** — `1` (names file) or `2` (interactive ranges).
   - Mode 1 then asks for the names file path (empty = `CONF/names.dat`).
   - Mode 2 then asks for ranges, one per line, e.g.:
     ```
     1 50 graphite
     51 80 bulk_hbn
     0 0 end
     ```

Example session:

```
$ ./make_ab_outcars
 Number of OUTCAR files found =  80
 CONF/OUTCAR1: atom types =    1, atoms =    64
    C  64
 ...
 Enter the number of extra structures to take back from the last one
 (0 = last structure only, n = last structure + n steps back):
0
 Choose how the system names are assigned:
   1 = one name per OUTCAR read from a names file
   2 = names given by ranges of OUTCAR numbers
2
 Enter the first and last OUTCAR number of a range and its name,
 e.g. "1 100 graphite" (the name must not contain blanks).
 Ranges can be entered repeatedly; enter "0 0 end" when done:
1 50 graphite
51 80 bulk_hbn
0 0 end
 ...
 Configurations found in the OUTCAR files =  80
 Configurations written to ML_AB          =  80
 The number of configurations in ML_AB (line 5) is already set
```

## Output

- **`ML_AB`** — the training file for `ML_MODE = select`, containing:
  - the correct configuration count (line 5)
  - the union of all atom types across every OUTCAR, with masses
  - placeholder reference atomic energies (`EATOM = 0.0`) — edit by hand if
    you have better values
  - placeholder basis-set entries (`1 1`) per atom type
  - the configurations themselves (lattice, positions, forces, stress),
    grouped by system name
- **`energy.dat`** — one total energy (eV) per line, in the same order the
  configurations were written to `ML_AB`.

## Notes and caveats

- `CTIFOR` is hardcoded to `0.002` in the source (`CTIFOR` parameter) —
  change and recompile if you need a different value.
- Atom types are read from the `VRHFIN` line of each POTCAR block in the
  OUTCAR, falling back to `TITEL` if `VRHFIN` isn't found.
- If the same atom type shows a different mass in different OUTCARs, the
  program prints a warning but keeps the mass from the first OUTCAR it saw
  that type in — it does not stop.
- Energy is taken from the last `free  energy` line before each `POSITION`
  block (i.e. `TOTEN`).
- Fixed limits in the source: `MXT = 50` (max atom types per OUTCAR),
  `LNM = 40` (max characters in a system name). Increase and recompile if
  you hit either limit.
- If an OUTCAR is missing energy, positions, stress, or lattice vectors, the
  program prints a warning naming the file rather than failing silently.

# make_ab_merge_part

Merge several VASP `ML_ABN`/`ML_AB` files (from different directories or
runs) into a single `ML_AB_merge` file, ready for refitting — no need to run
`ML_MODE = select` again unless you want to.

## What it does

- For each input, reads `DIR/ML_ABN` (falling back to `DIR/ML_AB` if
  `ML_ABN` doesn't exist), or reads the given path directly if it's already
  a file rather than a directory.
- For every file, prints how many configurations it has and lets you choose
  which ones to keep (a "first last" range, or `0 0`/empty = all) — so you
  can merge any subset of any file, not just the last one.
- Rebuilds the merged header from scratch: the number of configurations and
  the number of basis sets per atom type are **counted**, never copied from
  the input headers, and the maxima (atoms per system, atoms per atom type)
  are taken over all files.
- Merges atom types **by name** across files — the files don't need the same
  species, the same number of species, or one being a subset of another.
- Verifies each input file: the declared configuration count is checked
  against what's actually present, and every basis set is checked to point
  to a configuration that exists.
- Compares atomic masses and reference atomic energies for a given species
  across files and warns if they differ (the value from the first file that
  introduces that species is the one kept in the output).
- Configurations are copied **verbatim** — any line length, any number of
  atoms per configuration is handled — and renumbered sequentially in the
  merged file, with basis sets renumbered to match.

## Requirements

- `gfortran` (no external libraries)

## Compilation

```bash
gfortran -o make_ab_merge_part make_ab_merge_part.f90
```

## Usage

```bash
./make_ab_merge_part [dir1 dir2 ...]
```

- Each argument is **either**:
  - a directory — `DIR/ML_ABN` is read if it exists, otherwise `DIR/ML_AB`, or
  - a path to an `ML_AB`/`ML_ABN` file itself.
- **Without arguments**, the two directories hardcoded in the source
  (`DIRDEF = 'data1', 'data2'`) are used — edit `DIRDEF` and recompile if you
  want different defaults.
- The merge order follows the order of the arguments.
- Output is always written to `ML_AB_merge` in the current directory,
  **overwriting any existing file of that name without asking** — rename or
  copy it elsewhere once you're happy with it.
- If you give a single file, `ML_AB_merge` ends up as a (possibly filtered)
  copy of it; the program prints a warning to remind you of this.

Example:

```bash
./make_ab_merge_part C_N-C_i_0024 C_N-C_i_0025 C_N-C_i_0046
```

Example session:

```
$ ./make_ab_merge_part run1 run2

Reading run1/ML_ABN
 configurations =       40, species =    2: C  B
    basis sets per species =            5           3

Reading run2/ML_ABN
 configurations =       25, species =    2: C  N
    basis sets per species =            3           2

Merged species list:    3: C B N
Note: run1/ML_ABN does not contain all the species
Note: run2/ML_ABN does not contain all the species

File run1/ML_ABN: configurations =       40
Enter the first and last configuration to keep ("0 0" or empty = all): 

   keeping  40  configurations

File run2/ML_ABN: configurations =       25
Enter the first and last configuration to keep ("0 0" or empty = all): 1 20
   keeping  20  configurations

The number of configurations            =  60
The numbers of basis sets per atom type =            5           3           2
Written to ML_AB_merge
```

## Output

`ML_AB_merge`, with:

- the rebuilt header (correct configuration count; the union of atom types
  by name; recomputed maxima for atoms per system and per atom type)
- reference atomic energies and masses per species (from whichever file
  first introduced that species)
- basis sets renumbered to match the merged configuration numbering, listed
  per species in file order
- the kept configurations, renumbered sequentially, copied verbatim from
  their source files

## Checks and warnings

Hard errors (the program stops):
- an input file's declared configuration count doesn't match what's
  actually in the file
- a basis set points to a configuration index outside `1..NCONF` of its own
  file
- configurations in an input file aren't numbered `1, 2, 3, ...`
  consecutively
- an input `ML_AB`/`ML_ABN` file is missing a required header block (number
  of configurations, number/list of atom types, max atoms per system/type,
  or the basis-set section)

Warnings only (the program continues):
- a species has a different mass or reference energy in different files
- the version line (first line of the file) differs between files
- a file's basis sets are all placeholders (one `"1 1"` per species) — this
  usually means that file hasn't been through `ML_MODE = select` yet
- reference atomic energy or atomic mass missing from a file's header
  (treated as `0`)
- a file doesn't contain all the species seen in the merged set

## Notes and caveats

- Species names are stored as `LEN=8` strings (`LNAM` parameter) — longer
  names are truncated. Edit and recompile if you need longer names.
- Lines are read into `LEN=1024` buffers (`LLIN` parameter) — plenty for the
  usual one-atom-per-line `ML_AB` format, but edit and recompile if you ever
  hit longer lines.
- The output file name (`ML_AB_merge`) is fixed in the source (`OUTFILE`
  parameter).

## Typical workflow

This program is the natural follow-up to `make_ab_outcars.f90`: build an
initial `ML_AB` from a batch of OUTCARs, run VASP with `ML_MODE = select`
in each configuration's folder (which writes out an `ML_ABN` containing the
locally selected reference configurations), then use this program to fold
all those per-folder `ML_ABN` files back into one training set — ready to
refit without needing to re-run `ML_MODE = select`.

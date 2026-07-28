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

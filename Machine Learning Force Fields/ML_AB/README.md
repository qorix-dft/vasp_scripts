# Machine Learning Force Fields: ML_AB

Tools for VASP machine-learned force field (MLFF) training data.

## rm_structures.f90

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

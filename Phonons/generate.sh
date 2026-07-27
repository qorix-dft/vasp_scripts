#!/bin/bash
# generate_run_scripts.sh
#
# Generates a set of batch-scheduler run scripts for pure ab initio VASP
# relaxations (no ML force field) covering all POSCAR-* displacements,
# split into groups so each group is submitted as its own job.
#
# This version is scheduler-agnostic: edit the "SUPERCOMPUTER / SCHEDULER
# SETTINGS" block below to match whichever machine you're using (SLURM,
# PBS/Torque, LSF, ...). It's pre-filled with the Pawsey/SLURM settings
# currently in use -- replace that block wholesale for a different machine.
#
# Usage:
#   1. Edit the SUPERCOMPUTER / SCHEDULER SETTINGS block, and
#      TOTAL_POSCARS / GROUP_SIZE, below.
#   2. Run this script once from the directory containing POSCAR-*,
#      KPOINTS, POTCAR, INCAR:
#        bash generate_run_scripts.sh
#   3. Choose whether to submit the jobs:
#        bash generate_run_scripts.sh              generate, then asks y/n
#        bash generate_run_scripts.sh --submit      generate + submit all, no prompt
#        bash generate_run_scripts.sh --no-submit   generate only, no prompt
#
#      (Manual submission always still works too: <SUBMIT_CMD> run_01.sh, etc.)

# ---- Submission behaviour ----
# Used only when no --submit/--no-submit flag is given on the command line.
# "ask" -> prompt y/n after generating the files | "yes" -> auto-submit | "no" -> never auto-submit
SUBMIT_DEFAULT="ask"

SUBMIT_MODE="$SUBMIT_DEFAULT"
case "$1" in
    --submit)    SUBMIT_MODE="yes" ;;
    --no-submit) SUBMIT_MODE="no"  ;;
    "")          ;;  # keep SUBMIT_DEFAULT
    *)
        echo "Unknown option: $1"
        echo "Usage: $0 [--submit | --no-submit]"
        exit 1
        ;;
esac

# ============================================================
#  SUPERCOMPUTER / SCHEDULER SETTINGS
#  Replace this whole block with the settings for whichever machine
#  you're currently using. Pre-filled below with Pawsey/SLURM.
# ============================================================

# Command used to submit a job script: sbatch (SLURM), qsub (PBS/Torque), bsub (LSF), ...
SUBMIT_CMD="sbatch"

# Job name prefix -- the group number is appended automatically (g01, g02, ...)
JOB_NAME_PREFIX="phabinit"

# Scheduler directives + environment/module setup, written exactly as they
# should appear at the top of each generated run file. Use the placeholder
# __JOBNAME__ wherever the per-group job name should go.
# Nothing in this block is evaluated now -- it's copied in verbatim, so
# runtime variables like $SLURM_SUBMIT_DIR are written as literal text and
# only get evaluated later, when the generated job actually runs.
JOB_HEADER=$(cat << 'HEADER_EOF'
#!/bin/bash -l
#SBATCH --job-name=__JOBNAME__
#SBATCH --account=
#SBATCH --partition=work
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=115G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=name@emai.com

module load hdf5/1.14.5-parallel-api-v112 netlib-scalapack/2.2.0 fftw/3.3.10

# VASP 6.4.3 + VTST
export VASP=/software/vasp_std

# OpenMP settings
export OMP_NUM_THREADS=1

# MPI / Slingshot settings
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_VERBOSE=1
export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)

ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"
HEADER_EOF
)

# The command used to actually launch VASP for a single displacement.
# Typically references the $VASP variable exported above, plus whatever
# parallel-launch syntax your scheduler uses (srun, mpirun, aprun, ibrun...).
# Same rule as above: nothing here is evaluated now, it's copied verbatim.
LAUNCH_CMD=$(cat << 'LAUNCH_EOF'
  srun -N "$SLURM_JOB_NUM_NODES" \
       -n "$SLURM_NTASKS" \
       -c "$OMP_NUM_THREADS" \
       -m block:block:block \
       "$VASP" > vasp.out
LAUNCH_EOF
)

# ============================================================
#  GROUPING PARAMETERS
# ============================================================
TOTAL_POSCARS=1080     # total number of POSCAR-* files to process
GROUP_SIZE=60           # number of displacements handled per run file/job
OUTDIR="."              # where the generated run_XX.sh files are written

# ============================================================
#  Internal logic -- shouldn't normally need editing below this line
# ============================================================

NUM_GROUPS=$(( (TOTAL_POSCARS + GROUP_SIZE - 1) / GROUP_SIZE ))

PAD_WIDTH=${#NUM_GROUPS}
if [ "$PAD_WIDTH" -lt 2 ]; then
    PAD_WIDTH=2
fi

# Per-displacement directory setup + POSCAR handling, common to any
# scheduler. __START__/__END__/__LAUNCH__ are filled in per group below.
LOOP_TEMPLATE=$(cat << 'LOOP_EOF'
for i in {__START__..__END__}; do
  if [ "$i" -lt 10 ]; then
    pref=00$i
  elif [ "$i" -lt 100 ]; then
    pref=0$i
  else
    pref=$i
  fi
  TDIR=disp-$pref
  mkdir -p "$TDIR"
  cd "$TDIR" || exit 1
  mv "../POSCAR-$pref" POSCAR
  cp ../KPOINTS ../POTCAR ../INCAR .
__LAUNCH__
  cd ../
done
LOOP_EOF
)

CREATED_FILES=()

for (( g=1; g<=NUM_GROUPS; g++ )); do
    START=$(( (g-1)*GROUP_SIZE + 1 ))
    END=$(( g*GROUP_SIZE ))
    if [ "$END" -gt "$TOTAL_POSCARS" ]; then
        END=$TOTAL_POSCARS
    fi

    gpad=$(printf "%0${PAD_WIDTH}d" "$g")
    JOBNAME="${JOB_NAME_PREFIX}_g${gpad}"
    RUNFILE="${OUTDIR}/run_${gpad}.sh"

    HEADER_INST="${JOB_HEADER//__JOBNAME__/$JOBNAME}"

    BODY_INST="${LOOP_TEMPLATE//__START__/$START}"
    BODY_INST="${BODY_INST//__END__/$END}"
    BODY_INST="${BODY_INST//__LAUNCH__/$LAUNCH_CMD}"

    {
        printf '%s\n' "$HEADER_INST"
        printf '\n'
        printf '# This run file (group %s) covers displacements %s to %s\n' "$gpad" "$START" "$END"
        printf '%s\n' "$BODY_INST"
    } > "$RUNFILE"

    chmod +x "$RUNFILE"
    CREATED_FILES+=("$RUNFILE")
    echo "Created $RUNFILE (displacements $START-$END)"
done

echo ""
echo "Done: generated $NUM_GROUPS run file(s) in ${OUTDIR}/"

# ---- Submit, or not, based on SUBMIT_MODE ----
if [ "$SUBMIT_MODE" = "ask" ]; then
    read -rp "Submit all $NUM_GROUPS run script(s) now with $SUBMIT_CMD? [y/N] " REPLY
    case "$REPLY" in
        [yY]|[yY][eE][sS]) SUBMIT_MODE="yes" ;;
        *)                  SUBMIT_MODE="no"  ;;
    esac
fi

if [ "$SUBMIT_MODE" = "yes" ]; then
    if ! command -v "$SUBMIT_CMD" >/dev/null 2>&1; then
        echo "$SUBMIT_CMD not found on this system -- skipping submission."
        echo "Submit manually once you're on your cluster, e.g.: $SUBMIT_CMD run_01.sh"
    else
        echo "Submitting $NUM_GROUPS job(s) with $SUBMIT_CMD..."
        for f in "${CREATED_FILES[@]}"; do
            "$SUBMIT_CMD" "$f"
        done
    fi
else
    echo "Not submitting. Submit manually with, e.g.: $SUBMIT_CMD run_${gpad}.sh"
fi

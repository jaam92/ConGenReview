#!/bin/bash
#SBATCH --job-name=dfe_sims
#SBATCH --output=/scratch2/jazlynmo/bmc/logfiles/dfe_%A_%a.out
#SBATCH --error=/scratch2/jazlynmo/bmc/logfiles/dfe_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH -p qcb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-200%100
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazlynmo@usc.edu

# One array task = one replicate. Tasks 1-100 are 10 haplotypes,
# 101-200 are 20 haplotypes. If you change N_REPS or SAMPLE_SIZES,
# set --array=1-(N_REPS x number of sample sizes).
# "%100" = at most 100 tasks running at once; adjust to your queue.
# Timing test: ~2 min (10 haplotypes) and ~12 min (20 haplotypes) on
# 1 CPU (extra CPUs were not used), so 45 min leaves a safety margin.
#
# Before the first submission:
#   mkdir -p /scratch2/jazlynmo/bmc/logfiles

set -e

N_REPS=100
SAMPLE_SIZES="10 20"
NU=0.05
T=0.01
N_STARTS=10

SCRIPT=/scratch2/jazlynmo/bmc/scripts/dfe.py
RESULTS_DIR=/scratch2/jazlynmo/bmc/dfe_results

# Map the array index to a sample size and replicate number
task=$SLURM_ARRAY_TASK_ID
size_index=$(( (task - 1) / N_REPS + 1 ))
replicate=$(( (task - 1) % N_REPS + 1 ))
ns=$(echo $SAMPLE_SIZES | cut -d' ' -f"$size_index")
seed=$(( ns * 1000000 + replicate ))
output="${RESULTS_DIR}/ns${ns}/rep_$(printf '%03d' "$replicate").tsv"

# Skip replicates that already finished (makes resubmission safe)
if [ -s "$output" ]; then
    echo "Already done: $output"
    exit 0
fi

# Load the conda environment with dadi and nlopt (CARC's recommended
# pattern: no conda module; the shell hook makes "conda activate" work
# in batch jobs).
module purge
eval "$(conda shell.bash hook)"
conda activate dadi   # or the full path if created with -p

# One BLAS thread per process. Otherwise OpenBLAS tries to start one
# thread per core on the node (e.g. 64) and can hit process limits.
# dfe.py parallelizes the selection caches with --cpus instead.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Task $task: ns=$ns replicate=$replicate seed=$seed ($(date))"
python "$SCRIPT" "$ns" "$NU" "$T" "$replicate" "$seed" "$output" \
    --n-starts "$N_STARTS" --cpus "$SLURM_CPUS_PER_TASK"
echo "Finished $(date)"
echo "Elapsed: $((SECONDS / 60)) min $((SECONDS % 60)) s"

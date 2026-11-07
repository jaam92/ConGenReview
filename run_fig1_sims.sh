#!/bin/sh
#SBATCH --job-name=fig1_sim_100rep
#SBATCH --output=/scratch2/jazlynmo/bmc/logfiles/fig1_sim_%j.out
#SBATCH --error=/scratch2/jazlynmo/bmc/logfiles/fig1_sim_%j.err
#SBATCH --time=06:00:00
#SBATCH -p qcb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazlynmo@usc.edu

# NOTE: SLURM does not create the logfiles directory. Before the
# first submission, run:
#   mkdir -p /scratch2/jazlynmo/bmc/logfiles

set -e

# Paths
SIM_SCRIPT=/scratch2/jazlynmo/bmc/scripts/fig1_sims_final.R
PLOT_SCRIPT=/scratch2/jazlynmo/bmc/scripts/fig1_plot_sims.R
OUTDIR=/scratch2/jazlynmo/bmc/fig1_simulations

mkdir -p "$OUTDIR"

# Load R. Adjust to the module names on your cluster
# (check with: module spider r)
module purge
module load ver/2506  gcc/14.3.0 r/4.5.2

# One BLAS thread per worker; parallelism comes from the 16 R workers
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Tells both R scripts where to write/read CSVs and PDFs
export FIG1_OUTDIR="$OUTDIR"

echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Output dir: $OUTDIR"
echo "Started:    $(date)"

# 1) Simulate and write CSVs (skips simulations already saved)
Rscript "$SIM_SCRIPT"

# 2) Draw all figures from the CSVs
Rscript "$PLOT_SCRIPT"

echo "Finished:   $(date)"

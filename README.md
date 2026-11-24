# Code for Genomics in Conservation

This repository contains code used for the simulation-based analyses and figures presented in **“Genomics in Conservation: Analyzing and Interpreting Data for Reliable Action.”** The analyses illustrate how study design and analytical choices can affect inference in conservation genomics.

## Overview

The code includes simulations and visualizations examining:

- the joint effects of sample size, marker number, and population differentiation (`FST`) on the detection of population structure;
- error in demographic inference across sample sizes; and
- the effects of sample size on inference of the distribution of fitness effects.

These simulations are intended to illustrate general analytical behavior under controlled conditions. They do not reproduce all sources of complexity found in empirical genomic datasets, such as uneven missingness, ascertainment bias, variable sequencing depth, or heterogeneous linkage disequilibrium.

## Repository contents

| File | Purpose |
| --- | --- |
| `Figure1_sims_final.R` | Simulates population-genetic data across combinations of sample size, marker number, and `FST`, and calculates the power to detect population structure. |
| `run_fig1_sims.sh` | Shell wrapper for running the Figure 1 simulations. |
| `Figure1_detection_power.csv` | Detection-power results used to generate Figure 1. |
| `Figure1_plot_sims.R` | Summarizes and plots the Figure 1 simulation results. |
| `suppFig_LD_detection_power.csv` | Detection-power results for the supplementary analysis incorporating linked markers. |
| `Figure2_run_gone2.sh` | Runs the GONE2 demographic-inference analysis used for Figure 2. |
| `Figure2_run_msmc2.sh` | Runs the MSMC2 demographic-inference analysis used for Figure 2. |
| `Figure2_plot_msmc2_gone2.r` | Processes and plots the GONE2 and MSMC2 results for Figure 2. |
| `dfe.py` | Simulates and re-infers a one-population demographic model and gamma DFE. Each invocation runs one replicate and writes a two-line TSV file containing the header and results, including failed-fit diagnostics. |
| `run_dfe_array.sh` | SLURM array-job wrapper that runs 100 replicates for each of two sample sizes: 10 and 20 haplotypes (5 and 10 diploid individuals). |
| `dfe_all_replicates.csv` | Consolidated results from all 200 Figure 3 simulation replicates. |
| `fig3_plot_dfe.R` | Summarizes the replicate-level results and generates the three-panel Figure 3 showing inferred demographic histories, inferred DFE distributions, and the proportions of mutations in four selection classes. |

## Software requirements

The repository contains R, Python, and shell scripts. The R plotting and simulation code uses packages including:

```r
install.packages(c(
  "tidyverse",
  "scales",
  "patchwork"
))
```

The Figure 3 workflow requires Python 3 with `dadi`, `nlopt`, and `numpy`. Figure 2 additionally requires GONE2 and MSMC2. Paths to external executables, scripts, results directories, and input files may need to be adjusted for the user's computing environment.

## Running the analyses

1. Clone or download this repository.
2. Open R or RStudio and set the repository as the working directory.
3. Install the required packages.
4. Run the scripts for the desired figure in the order shown below.

### Figure 1

Run the simulations and then generate the plot:

```bash
sbatch run_fig1_sims.sh
Rscript Figure1_plot_sims.R
```

If `Figure1_detection_power.csv` has already been generated, only the plotting script is needed:

```bash
Rscript Figure1_plot_sims.R
```

### Figure 2

Run GONE2 and MSMC2, then process and plot their outputs:

```bash
sbatch Figure2_run_gone2.sh
sbatch Figure2_run_msmc2.sh
Rscript Figure2_plot_msmc2_gone2.r
```

### Figure 3

Before submitting the simulations, edit `SCRIPT`, `RESULTS_DIR`, the SLURM log paths, partition, and email address in `run_dfe_array.sh` for the local computing environment. The default array contains 200 tasks: tasks 1–100 analyze 10 haplotypes and tasks 101–200 analyze 20 haplotypes.

```bash
sbatch run_dfe_array.sh
```

Each array task invokes `dfe.py` once and writes one replicate-level TSV file under `RESULTS_DIR/ns10/` or `RESULTS_DIR/ns20/`. Completed files are skipped when the array is resubmitted. To summarize the replicate files and generate Figure 3, set the input and output directories and run:

```bash
DFE_RESULTS_DIR=/path/to/dfe_results \
DFE_PLOT_DIR=/path/to/figure_outputs \
Rscript fig3_plot_dfe.R
```

The plotting script writes:

- `Fig3_demography_DFE.pdf`, the final three-panel Figure 3;
- `fig3_run_summary.csv`, a summary of successful fits, failed fits, and estimates at parameter bounds; and
- `fig3_error_summary.csv`, RMAE and plain-scale summaries of the inferred demographic and DFE parameters.

The repository also includes `dfe_all_replicates.csv`, which consolidates the 200 replicate-level results into a single file for inspection and downstream analysis.

The shell scripts may contain paths or computing-cluster settings that must be modified before use. The scripts set random-number seeds where appropriate to support reproducibility. Because some analyses include many simulation replicates, large numbers of markers, or computationally intensive demographic inference, run times and memory requirements will vary across computers. We use dadi version 2.4.5.

## Population-structure simulations

The population-structure simulations generate diploid genotypes for two populations while varying:

- the number of independent SNPs;
- `FST`; and
- the number of sampled haplotypes per population.

The two sample-size conditions correspond to 20 haplotypes per population (10 diploid individuals) and 8 haplotypes per population (4 diploid individuals). Principal component analysis is used to evaluate population separation across parameter combinations. Replicated simulations summarize the relationship among marker density, sample size, and population differentiation rather than relying on a single stochastic realization.

The primary simulations assume independent markers. In empirical datasets, linkage disequilibrium, missingness, ascertainment bias, and limited marker informativeness can reduce the effective amount of independent information available for detecting population structure. Results from the supplementary linked-marker analysis are provided in `suppFig_LD_detection_power.csv`.

## Demographic and DFE analyses

Figure 2 compares demographic histories inferred using GONE2 and MSMC2. The shell scripts run the respective programs, and `Figure2_plot_msmc2_gone2.r` processes their outputs for visualization.

Figure 3 examines how sample size affects demographic and distribution-of-fitness-effects inference. The simulations use an ancestral effective population size of 1,000 followed by a recent contraction to 5% of that size. The true DFE is gamma distributed with shape parameter alpha = 0.2 and scale parameter beta = 1,000.

For each sample size, `dfe.py` simulates a folded synonymous SFS under the true demography, infers demographic parameters from that SFS, and estimates the synonymous mutation-rate scaling parameter. It then simulates a folded nonsynonymous SFS under the true demography and infers the DFE using the fitted demography and estimated nonsynonymous scaling parameter. Both demographic and DFE optimization use multiple random starting points; failed fits and estimates near parameter bounds are retained and flagged rather than silently discarded.

The simulations compare 10 and 20 haplotypes, corresponding to 5 and 10 diploid individuals. `fig3_plot_dfe.R` displays replicate-level demographic histories and summarizes the median and 95% interval of the inferred DFE. It also translates the fitted DFE into the proportions of mutations that are effectively neutral, weakly deleterious, moderately deleterious, or strongly deleterious according to 2Nₐ|s|.

## Reproducibility and interpretation

Simulation settings and random-number seeds are specified directly in the scripts. Users can modify these values to explore alternative sample sizes, marker counts, levels of population differentiation, or numbers of simulation replicates.

The code is provided to support transparency and reproducibility. Results should be interpreted in relation to the assumptions of each simulation rather than as universal thresholds for conservation-genomic studies.

## Use of generative AI

Generative AI tools, including OpenAI ChatGPT, were used to assist with portions of code generation, revision, documentation, and debugging. All AI-assisted code was reviewed, edited, and validated by the authors before use. The authors take full responsibility for the implementation, interpretation, and reported results.

## Repository

<https://github.com/jaam92/ConGenReview>

## Questions

Questions about the analyses or code can be submitted through the repository's issue tracker.

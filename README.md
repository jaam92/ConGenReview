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
| `Figures3And4_dfe.py` | Generates or processes the demographic and distribution-of-fitness-effects analyses used for Figures 3 and 4. |
| `Figures3And4_run_dfe.sh` | Shell wrapper for running the Figures 3 and 4 analyses. |
| `Figures3And4_plot_error.R` | Calculates, summarizes, and plots inference error for Figures 3 and 4. |

## Software requirements

The repository contains R, Python, and shell scripts. The R plotting and simulation code uses packages including:

```r
install.packages(c(
  "tidyverse",
  "scales"
))
```

The analyses for Figures 2–4 also require the external software and Python packages referenced in the corresponding scripts, including GONE2 and MSMC2 for Figure 2. Paths to external executables and input files may need to be adjusted for the user's computing environment.

## Running the analyses

1. Clone or download this repository.
2. Open R or RStudio and set the repository as the working directory.
3. Install the required packages.
4. Run the scripts for the desired figure in the order shown below.

### Figure 1

Run the simulations and then generate the plot:

```bash
bash run_fig1_sims.sh
Rscript Figure1_plot_sims.R
```

If `Figure1_detection_power.csv` has already been generated, only the plotting script is needed:

```bash
Rscript Figure1_plot_sims.R
```

### Figure 2

Run GONE2 and MSMC2, then process and plot their outputs:

```bash
bash Figure2_run_gone2.sh
bash Figure2_run_msmc2.sh
Rscript Figure2_plot_msmc2_gone2.r
```

### Figures 3 and 4

Run the DFE workflow and then plot inference error:

```bash
bash Figures3And4_run_dfe.sh
Rscript Figures3And4_plot_error.R
```

The shell scripts may contain paths or computing-cluster settings that must be modified before use. The scripts set random-number seeds where appropriate to support reproducibility. Because some analyses include many simulation replicates, large numbers of markers, or computationally intensive demographic inference, run times and memory requirements will vary across computers.

## Population-structure simulations

The population-structure simulations generate diploid genotypes for two populations while varying:

- the number of independent SNPs;
- `FST`; and
- the number of sampled haplotypes per population.

The two sample-size conditions correspond to 20 haplotypes per population (10 diploid individuals) and 8 haplotypes per population (4 diploid individuals). Principal component analysis is used to evaluate population separation across parameter combinations. Replicated simulations summarize the relationship among marker density, sample size, and population differentiation rather than relying on a single stochastic realization.

The primary simulations assume independent markers. In empirical datasets, linkage disequilibrium, missingness, ascertainment bias, and limited marker informativeness can reduce the effective amount of independent information available for detecting population structure. Results from the supplementary linked-marker analysis are provided in `suppFig_LD_detection_power.csv`.

## Demographic and DFE analyses

Figure 2 compares demographic histories inferred using GONE2 and MSMC2. The shell scripts run the respective programs, and `Figure2_plot_msmc2_gone2.r` processes their outputs for visualization.

Figures 3 and 4 examine how sample size affects demographic and distribution-of-fitness-effects inference. `Figures3And4_dfe.py` and `Figures3And4_run_dfe.sh` implement the analysis workflow, while `Figures3And4_plot_error.R` summarizes deviations from the simulated values and generates the corresponding figures.

## Reproducibility and interpretation

Simulation settings and random-number seeds are specified directly in the scripts. Users can modify these values to explore alternative sample sizes, marker counts, levels of population differentiation, or numbers of simulation replicates.

The code is provided to support transparency and reproducibility. Results should be interpreted in relation to the assumptions of each simulation rather than as universal thresholds for conservation-genomic studies.

## Use of generative AI

Generative AI tools, including OpenAI ChatGPT, were used to assist with portions of code generation, revision, documentation, and debugging. All AI-assisted code was reviewed, edited, and validated by the authors before use. The authors take full responsibility for the implementation, interpretation, and reported results.

## Repository

<https://github.com/jaam92/ConGenReview>

## Questions

Questions about the analyses or code can be submitted through the repository's issue tracker.

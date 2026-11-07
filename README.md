# Code for Genomics in Conservation

This repository contains code used for the simulation-based analyses and figures presented in **“Genomics in Conservation: Analyzing and Interpreting Data for Reliable Action.”** The analyses illustrate how study design and analytical choices can affect inference in conservation genomics.

## Overview

The code includes simulations and visualizations examining:

- the joint effects of sample size, marker number, and population differentiation (`FST`) on the detection of population structure;
- error in demographic inference across sample sizes; and
- the effects of sample size on inference of the distribution of fitness effects.

These simulations are intended to illustrate general analytical behavior under controlled conditions. They do not reproduce all sources of complexity found in empirical genomic datasets, such as uneven missingness, ascertainment bias, variable sequencing depth, or heterogeneous linkage disequilibrium.

## Software requirements

The analyses were conducted in R. Required packages include:

```r
install.packages(c(
  "tidyverse",
  "scales"
))
```

Additional package requirements, when applicable, are listed within the individual scripts.

## Running the analyses

1. Clone or download this repository.
2. Open R or RStudio and set the repository as the working directory.
3. Install the required packages.
4. Run each analysis script independently.

For example:

```r
source("script_name.R")
```

The scripts set random-number seeds where appropriate to support reproducibility. Because some analyses include many simulation replicates or large numbers of markers, run times and memory requirements may vary across computers.

## Population-structure simulations

The population-structure simulations generate diploid genotypes for two populations while varying:

- the number of independent SNPs;
- `FST`; and
- the number of sampled haplotypes per population.

The two sample-size conditions correspond to 20 haplotypes per population (10 diploid individuals) and 8 haplotypes per population (4 diploid individuals). Principal component analysis is used to evaluate population separation across parameter combinations. Replicated simulations summarize the relationship among marker density, sample size, and population differentiation rather than relying on a single stochastic realization.

The simulations assume independent markers. In empirical datasets, linkage disequilibrium, missingness, ascertainment bias, and limited marker informativeness can reduce the effective amount of independent information available for detecting population structure.

## Reproducibility and interpretation

Simulation settings and random-number seeds are specified directly in the scripts. Users can modify these values to explore alternative sample sizes, marker counts, levels of population differentiation, or numbers of simulation replicates.

The code is provided to support transparency and reproducibility. Results should be interpreted in relation to the assumptions of each simulation rather than as universal thresholds for conservation-genomic studies.

## Use of generative AI

Generative AI tools, including OpenAI ChatGPT, were used to assist with portions of code generation, revision, documentation, and debugging. All AI-assisted code was reviewed, edited, and validated by the authors. The authors take full responsibility for the implementation, interpretation, and reported results.

## Repository

<https://github.com/jaam92/ConGenReview>

## Questions

Questions about the analyses or code can be submitted through the repository's issue tracker.

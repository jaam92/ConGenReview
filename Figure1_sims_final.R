library(tidyverse)

#==========================================================
# Fig. 1 simulations: run simulations and write CSV files
#
# Outputs (all in out_dir):
#   raw_main_results.csv, raw_ld_results.csv
#       every replicate; also used as a cache, so a rerun loads
#       these instead of re-simulating
#   one summary CSV per plot (read by fig1_plot.R):
#       fig_correlation.csv          |r| vs SNPs (+ null line)
#       fig_detection_power.csv      detection-power heatmap
#       fig_collapse.csv             collapse plot
#       suppFig_separation.csv          standardized separation (+ null)
#       suppFig_PVE1.csv                variance explained (+ null)
#       suppFig_LD_correlation.csv      LD |r| line plot
#       suppFig_LD_detection_power.csv  LD detection-power heatmap
#       suppFig_example_PCA_scores.csv  PC1 / PC2 example plots
#==========================================================


#----------------------------------------------------------
# Settings
#----------------------------------------------------------

use_parallel <- TRUE

# Workers = CPUs actually available (respects SLURM batch and
# salloc allocations, cgroups, or local cores)
n_workers <- max(1, parallelly::availableCores())

# Output directory: set by the SLURM wrapper, else current dir
out_dir <- Sys.getenv("FIG1_OUTDIR", unset = getwd())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

# Raw results are cached as CSV. Set TRUE to force re-simulating.
force_rerun <- FALSE

# One BLAS/OpenMP thread per worker. Env vars affect workers
# started later; RhpcBLASctl (if installed) also sets this session.
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

message("Workers: ", n_workers, " | Output directory: ", out_dir)

# Replicates. Detection rates near the threshold are sensitive
# to the null cutoff (95th percentile of |r| with FST = 0), so
# the null needs many more runs than the structured conditions.
# Conditions with few (effective) SNPs are cheap and are where
# detection is uncertain, so they get more replicates; large-SNP
# conditions are expensive and sit at 100%, so they get fewer.
# For a quick timing test, divide all four by 10.
n_reps_small     <- 500   # structured runs, <= 1,000 (effective) SNPs
n_reps_large     <- 100   # structured runs, > 1,000 (effective) SNPs
n_reps_null_small <- 2000 # null runs, <= 1,000 (effective) SNPs
n_reps_null_large <- 500  # null runs, > 1,000 (effective) SNPs

reps_for <- function(fst, n_loci) {
  small <- n_loci <= 1000
  ifelse(fst == 0,
         ifelse(small, n_reps_null_small, n_reps_null_large),
         ifelse(small, n_reps_small, n_reps_large))
}


#----------------------------------------------------------
# Parameters (used for every comparison)
#----------------------------------------------------------

FSTs <- c(0.025, 0.050, 0.100, 0.125, 0.150, 0.200)
n_haplotypes <- c(4, 8, 12, 16, 20, 40)
Ls <- c(10, 50, 100, 1000, 10000, 100000)

# LD simulation
ld_L_totals <- c(1000, 10000, 100000)
ld_block_sizes <- c(1, 10, 100)
ld_r2 <- c(1, 0.5, 0.1)

# Example PCAs
example_Ls <- c(10, 100, 1000, 10000)


#==========================================================
# Core helpers
#==========================================================

chunk_sizes <- function(total, size) {
  sizes <- rep(size, total %/% size)
  if (total %% size > 0) sizes <- c(sizes, total %% size)
  sizes
}

draw_freqs <- function(n_loci, fst) {
  p0 <- runif(n_loci, 0.1, 0.9)
  delta <- sqrt(fst * p0 * (1 - p0))
  list(
    p1 = pmin(pmax(p0 + delta, 0), 1),
    p2 = pmin(pmax(p0 - delta, 0), 1)
  )
}

# Build ZZ' chunk by chunk (Z = centered, scaled genotypes,
# invariant SNPs removed). PC1 scores and PVE are identical to
# prcomp(center = TRUE, scale. = TRUE).
accumulate_K <- function(gen_chunk, sizes, n) {
  K <- matrix(0, n, n)
  n_var <- 0L
  for (cs in sizes) {
    X <- gen_chunk(cs)
    m <- colMeans(X)
    s <- sqrt(colSums((X - rep(m, each = n))^2) / (n - 1))
    keep <- s > 0
    if (any(keep)) {
      Z <- (X[, keep, drop = FALSE] - rep(m[keep], each = n)) /
        rep(s[keep], each = n)
      K <- K + tcrossprod(Z)
      n_var <- n_var + sum(keep)
    }
  }
  list(K = K, n_var = n_var)
}

pc_scores <- function(K, n_pcs = 1) {
  e <- eigen(K, symmetric = TRUE)
  ev <- pmax(e$values, 0)
  scores <- sweep(
    e$vectors[, seq_len(n_pcs), drop = FALSE], 2,
    sqrt(ev[seq_len(n_pcs)]), "*"
  )
  list(scores = scores, pve = ev / sum(ev))
}

pca_stats <- function(acc, pop) {
  if (acc$n_var < 2) {
    return(tibble(
      n_variable_snps = acc$n_var,
      PVE1 = NA_real_,
      pc1_separation = NA_real_,
      pc1_pop_correlation = NA_real_
    ))
  }

  pcs <- pc_scores(acc$K, 1)
  pc1 <- pcs$scores[, 1]

  g1 <- pc1[pop == 0]
  g2 <- pc1[pop == 1]
  n1 <- length(g1)
  n2 <- length(g2)
  pooled_sd <- sqrt(((n1 - 1) * var(g1) + (n2 - 1) * var(g2)) / (n1 + n2 - 2))

  tibble(
    n_variable_snps = acc$n_var,
    PVE1 = pcs$pve[1],
    pc1_separation = if (is.finite(pooled_sd) && pooled_sd > 0) {
      abs(mean(g1) - mean(g2)) / pooled_sd
    } else {
      NA_real_
    },
    pc1_pop_correlation = abs(cor(pc1, pop))
  )
}

# Parallel runner. Biggest jobs first + many small chunks keeps
# all workers busy until the end.
run_grid <- function(grid, fun, seed) {
  if (use_parallel && requireNamespace("furrr", quietly = TRUE)) {
    # Batch (Rscript): forked workers. Interactive: separate R
    # sessions, because forking an interactive session whose
    # BLAS/OpenMP threads are running can crash workers.
    if (!interactive() && future::supportsMulticore()) {
      future::plan(future::multicore, workers = n_workers)
    } else {
      future::plan(future::multisession, workers = n_workers)
    }
    on.exit(future::plan(future::sequential), add = TRUE)
    res <- furrr::future_pmap(
      grid, fun,
      .options = furrr::furrr_options(seed = seed, scheduling = 10)
    )
  } else {
    if (use_parallel) message("Install 'furrr' to run in parallel; running sequentially.")
    set.seed(seed)
    res <- pmap(grid, fun, .progress = TRUE)
  }
  bind_rows(res)
}

# Load cached raw results (CSV) if present; otherwise run and save
run_cached <- function(cache_file, grid, fun, seed) {
  if (!force_rerun && file.exists(cache_file)) {
    message("Loading cached results from ", cache_file)
    return(read_csv(cache_file, show_col_types = FALSE))
  }
  t0 <- Sys.time()
  message("Running ", nrow(grid), " simulations -> ", cache_file)
  res <- run_grid(grid, fun, seed)
  write_csv(res, cache_file)
  message("Finished in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
  res
}


#==========================================================
# Genotype generators
#==========================================================

gen_independent <- function(n_dip, fst) {
  function(n_snps) {
    f <- draw_freqs(n_snps, fst)
    N <- n_dip * n_snps
    P1 <- rep(f$p1, each = n_dip)
    P2 <- rep(f$p2, each = n_dip)
    rbind(
      matrix((runif(N) < P1) + (runif(N) < P1), nrow = n_dip),
      matrix((runif(N) < P2) + (runif(N) < P2), nrow = n_dip)
    )
  }
}

# Linked SNPs: each SNP on each haplotype copies its block's
# founder allele with probability c, otherwise it is redrawn at
# the same frequency. FST is unchanged, and the expected
# within-population r^2 between SNPs in a block is c^4.
gen_linked <- function(n_dip, fst, block_size, copy_prob) {

  sim_pop <- function(p) {
    n_hap <- 2 * n_dip
    n_blocks <- length(p)
    n_snps <- n_blocks * block_size

    founder <- matrix(runif(n_hap * n_blocks) < rep(p, each = n_hap), nrow = n_hap)
    H <- matrix(
      runif(n_hap * n_snps) < rep(rep(p, each = block_size), each = n_hap),
      nrow = n_hap
    )
    copy <- matrix(runif(n_hap * n_snps) < copy_prob, nrow = n_hap)
    founder_exp <- founder[, rep(seq_len(n_blocks), each = block_size), drop = FALSE]
    H[copy] <- founder_exp[copy]

    H[seq(1, n_hap, 2), , drop = FALSE] + H[seq(2, n_hap, 2), , drop = FALSE]
  }

  function(n_blocks) {
    f <- draw_freqs(n_blocks, fst)
    rbind(sim_pop(f$p1), sim_pop(f$p2))
  }
}


#==========================================================
# Analysis wrappers
#==========================================================

analyze_independent <- function(fst, L, n_haplotypes_per_pop, replicate_id) {
  n_dip <- n_haplotypes_per_pop / 2
  n <- 2 * n_dip
  sizes <- chunk_sizes(L, max(1, floor(4e6 / n)))
  acc <- accumulate_K(gen_independent(n_dip, fst), sizes, n)

  bind_cols(
    tibble(fst = fst, L = L, n_haplotypes_per_pop = n_haplotypes_per_pop,
           replicate_id = replicate_id),
    pca_stats(acc, rep(c(0, 1), each = n_dip))
  )
}

analyze_linked <- function(fst, L_total, block_size, r2,
                           n_haplotypes_per_pop, replicate_id) {
  n_dip <- n_haplotypes_per_pop / 2
  n <- 2 * n_dip
  n_blocks <- L_total %/% block_size
  sizes <- chunk_sizes(n_blocks, max(1, floor(1e6 / (n * block_size))))
  acc <- accumulate_K(gen_linked(n_dip, fst, block_size, r2^(1 / 4)), sizes, n)

  bind_cols(
    tibble(fst = fst, L_total = L_total, block_size = block_size, r2 = r2,
           n_haplotypes_per_pop = n_haplotypes_per_pop, replicate_id = replicate_id),
    pca_stats(acc, rep(c(0, 1), each = n_dip))
  )
}

example_scores <- function(fst, L, n_haplotypes_per_pop) {
  n_dip <- n_haplotypes_per_pop / 2
  n <- 2 * n_dip
  acc <- accumulate_K(gen_independent(n_dip, fst),
                      chunk_sizes(L, max(1, floor(4e6 / n))), n)
  pcs <- pc_scores(acc$K, 2)
  tibble(
    fst = fst, L = L, n_haplotypes_per_pop = n_haplotypes_per_pop,
    individual = seq_len(n),
    PC1 = pcs$scores[, 1], PC2 = pcs$scores[, 2],
    pop = rep(c("population 1", "population 2"), each = n_dip)
  )
}


#==========================================================
# Part A: independent SNPs (+ FST = 0 null)
#==========================================================

main_grid <- crossing(
  fst = c(0, FSTs),
  L = Ls,
  n_haplotypes_per_pop = n_haplotypes
) %>%
  mutate(n_rep = reps_for(fst, L)) %>%
  uncount(n_rep, .id = "replicate_id") %>%
  arrange(desc(L * n_haplotypes_per_pop))

simulation_results <- run_cached("raw_main_results.csv", main_grid,
                                 analyze_independent, seed = 120)

null_summary <- simulation_results %>%
  filter(fst == 0) %>%
  group_by(L, n_haplotypes_per_pop) %>%
  summarize(
    n_null_replicates = n(),
    null_median      = median(pc1_pop_correlation, na.rm = TRUE),
    null_q95         = quantile(pc1_pop_correlation, 0.95, na.rm = TRUE),
    null_sep_q95     = quantile(pc1_separation, 0.95, na.rm = TRUE),
    null_PVE1_median = median(PVE1, na.rm = TRUE),
    .groups = "drop"
  )

structured_results <- simulation_results %>% filter(fst > 0)

simulation_summary <- structured_results %>%
  group_by(fst, L, n_haplotypes_per_pop) %>%
  summarize(
    median_correlation = median(pc1_pop_correlation, na.rm = TRUE),
    correlation_lower  = quantile(pc1_pop_correlation, 0.025, na.rm = TRUE),
    correlation_upper  = quantile(pc1_pop_correlation, 0.975, na.rm = TRUE),
    median_separation  = median(pc1_separation, na.rm = TRUE),
    separation_lower   = quantile(pc1_separation, 0.025, na.rm = TRUE),
    separation_upper   = quantile(pc1_separation, 0.975, na.rm = TRUE),
    median_PVE1        = median(PVE1, na.rm = TRUE),
    PVE1_lower         = quantile(PVE1, 0.025, na.rm = TRUE),
    PVE1_upper         = quantile(PVE1, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(null_summary, by = c("L", "n_haplotypes_per_pop")) %>%
  mutate(
    # total individuals = 2 pops x (haplotypes / 2) = n_haplotypes_per_pop
    bbp_scale = fst * sqrt(n_haplotypes_per_pop * L),
    adjusted_correlation = (median_correlation - null_median) / (1 - null_median)
  )

power_summary <- structured_results %>%
  left_join(select(null_summary, L, n_haplotypes_per_pop, null_q95, n_null_replicates),
            by = c("L", "n_haplotypes_per_pop")) %>%
  group_by(fst, L, n_haplotypes_per_pop) %>%
  summarize(
    power = mean(coalesce(pc1_pop_correlation > null_q95, FALSE)),
    n_replicates = n(),
    null_q95 = first(null_q95),
    n_null_replicates = first(n_null_replicates),
    .groups = "drop"
  ) %>%
  mutate(predicted_detectable = fst * sqrt(n_haplotypes_per_pop * L) > 1)


#==========================================================
# Part B: linked SNPs (+ FST = 0 null for each LD design)
#==========================================================

ld_grid <- crossing(
  fst = c(0, FSTs),
  L_total = ld_L_totals,
  block_size = ld_block_sizes,
  r2 = ld_r2,
  n_haplotypes_per_pop = n_haplotypes
) %>%
  filter(block_size > 1 | r2 == 1) %>%  # r^2 meaningless for 1 SNP per block
  mutate(n_rep = reps_for(fst, L_total / block_size)) %>%  # by number of blocks
  uncount(n_rep, .id = "replicate_id") %>%
  arrange(desc(L_total * n_haplotypes_per_pop))

ld_results <- run_cached("raw_ld_results.csv", ld_grid, analyze_linked, seed = 121)

# The null depends on the LD design (e.g. 1,000 SNPs in blocks of
# 100 behave like 10 loci), so each design gets its own null
ld_null <- ld_results %>%
  filter(fst == 0) %>%
  group_by(L_total, block_size, r2, n_haplotypes_per_pop) %>%
  summarize(null_q95 = quantile(pc1_pop_correlation, 0.95, na.rm = TRUE),
            n_null_replicates = n(),
            .groups = "drop")

ld_summary <- ld_results %>%
  filter(fst > 0) %>%
  group_by(fst, L_total, block_size, r2, n_haplotypes_per_pop) %>%
  summarize(
    median_correlation = median(pc1_pop_correlation, na.rm = TRUE),
    correlation_lower  = quantile(pc1_pop_correlation, 0.025, na.rm = TRUE),
    correlation_upper  = quantile(pc1_pop_correlation, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(n_blocks = L_total / block_size)

ld_power_summary <- ld_results %>%
  filter(fst > 0) %>%
  left_join(ld_null, by = c("L_total", "block_size", "r2", "n_haplotypes_per_pop")) %>%
  group_by(fst, L_total, block_size, r2, n_haplotypes_per_pop) %>%
  summarize(
    power = mean(coalesce(pc1_pop_correlation > null_q95, FALSE)),
    n_replicates = n(),
    null_q95 = first(null_q95),
    n_null_replicates = first(n_null_replicates),
    .groups = "drop"
  ) %>%
  mutate(
    n_blocks = L_total / block_size,
    # Theory, counting each linkage block as one independent locus
    predicted_detectable = fst * sqrt(n_haplotypes_per_pop * n_blocks) > 1
  )


#==========================================================
# Part C: example PCAs
#==========================================================

set.seed(122)
example_pca <- crossing(
  fst = c(0, FSTs),
  L = example_Ls,
  n_haplotypes_per_pop = n_haplotypes
) %>%
  pmap(example_scores) %>%
  bind_rows()


#==========================================================
# Write one summary CSV per plot
#==========================================================

write_csv(
  select(simulation_summary, fst, L, n_haplotypes_per_pop,
         median_correlation, correlation_lower, correlation_upper, null_q95),
  "fig_correlation.csv"
)
write_csv(power_summary, "fig_detection_power.csv")
write_csv(
  select(simulation_summary, fst, L, n_haplotypes_per_pop,
         bbp_scale, median_correlation, null_median, adjusted_correlation),
  "fig_collapse.csv"
)
write_csv(
  select(simulation_summary, fst, L, n_haplotypes_per_pop,
         median_separation, separation_lower, separation_upper, null_sep_q95),
  "suppFig_separation.csv"
)
write_csv(
  select(simulation_summary, fst, L, n_haplotypes_per_pop,
         median_PVE1, PVE1_lower, PVE1_upper, null_PVE1_median),
  "suppFig_PVE1.csv"
)
write_csv(ld_summary, "suppFig_LD_correlation.csv")
write_csv(ld_power_summary, "suppFig_LD_detection_power.csv")
write_csv(example_pca, "suppFig_example_PCA_scores.csv")

message("Summary CSVs written to ", out_dir)

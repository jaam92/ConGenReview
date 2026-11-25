library(tidyverse)
library(patchwork)

#==========================================================
# Figure 3: demographic and DFE inference at two sample sizes
#
#  A. Inferred population-size histories (one facet per sample size)
#  B. Inferred DFE: density per log10 |s| (median + 95% interval)
#  C. Share of new mutations in selection classes (2Na|s|)
#
# Reads the per-replicate TSVs written by dfe.py
# (results_dir/ns*/rep_*.tsv). Also writes CSV summaries.
#==========================================================

# Folder containing ns10/ and ns20/ (searched recursively for rep_*.tsv).
# Edit this path, or set DFE_RESULTS_DIR before running.
results_dir <- Sys.getenv("DFE_RESULTS_DIR",
                          unset = "/scratch2/jazlynmo/bmc/dfe_results")
out_dir <- Sys.getenv("DFE_PLOT_DIR", unset = results_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!interactive()) pdf(NULL)

N_anc <- 1000  # ancestral Ne used in the simulations


#----------------------------------------------------------
# Read results
#----------------------------------------------------------

# If the combined file already exists, read it (fast, and it can be
# copied to a laptop on its own). Otherwise read the per-replicate
# TSVs and create it. Set rebuild_combined <- TRUE after adding or
# rerunning replicates so the combined file is regenerated.
rebuild_combined <- FALSE
combined_file <- file.path(out_dir, "dfe_all_replicates.csv")

result_col_types <- cols(status = col_character(), error_message = col_character(),
                         .default = col_double())

if (!rebuild_combined && file.exists(combined_file)) {
  message("Reading combined results from ", combined_file)
  res <- read_csv(combined_file, col_types = result_col_types, na = c("", "NA"))
} else {
  files <- list.files(results_dir, pattern = "^rep_.*\\.tsv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    stop("No dfe_all_replicates.csv in ", normalizePath(out_dir, mustWork = FALSE),
         " and no rep_*.tsv files under ", normalizePath(results_dir, mustWork = FALSE),
         ". Check results_dir / DFE_PLOT_DIR at the top of the script.")
  }
  message("Reading ", length(files), " replicate files from ", results_dir)
  res <- map_dfr(files, ~ read_tsv(.x, col_types = result_col_types, na = c("", "NA")))
  write_csv(res, combined_file)
  message("Wrote combined results to ", combined_file)
}

message("Replicates loaded: ", nrow(res))

# True values come from the files, so they always match the run
nu_true <- first(res$nu_true)
T_true <- first(res$T_true)
alpha_true <- first(res$alpha_true)
beta_true <- first(res$beta_true)

sizes <- sort(unique(res$sample_haplotypes), decreasing = TRUE)
size_long <- paste0(sizes, " haplotypes (", sizes / 2, " diploid individuals)")
size_short <- paste0("n = ", sizes)

# Colors: largest sample blue, smallest orange (colorblind-safe)
size_colors <- setNames(c("#0072B2", "#D55E00")[seq_along(sizes)], size_short)

res <- res %>%
  mutate(
    size_long = factor(paste0(sample_haplotypes, " haplotypes (",
                              sample_haplotypes / 2, " diploid individuals)"),
                       levels = size_long),
    size_short = factor(paste0("n = ", sample_haplotypes), levels = size_short)
  )


#----------------------------------------------------------
# Run summary and error table (for text / caption / supplement)
#----------------------------------------------------------

run_summary <- res %>%
  group_by(sample_haplotypes) %>%
  summarize(
    replicates_run = n(),
    converged = sum(status == "ok"),
    failed = sum(status != "ok"),
    demography_at_bound = sum(demog_at_bound == 1, na.rm = TRUE),
    dfe_at_bound = sum(dfe_at_bound == 1, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(run_summary, file.path(out_dir, "fig3_run_summary.csv"))
print(run_summary)

ok <- filter(res, status == "ok")

rmae <- function(fit, truth) mean(abs(fit - truth) / truth)
error_table <- ok %>%
  group_by(sample_haplotypes) %>%
  summarize(
    alpha_RMAE = rmae(alpha_fit, alpha_true),
    beta_RMAE = rmae(beta_fit, beta_true),
    contraction_size_RMAE = rmae(nu_fit, nu_true),
    contraction_time_RMAE = rmae(T_fit, T_true),
    # In plain units, for the text
    Ne_after_median = median(nu_fit * N_anc),
    Ne_after_2.5 = quantile(nu_fit * N_anc, 0.025),
    Ne_after_97.5 = quantile(nu_fit * N_anc, 0.975),
    generations_median = median(T_fit * 2 * N_anc),
    generations_2.5 = quantile(T_fit * 2 * N_anc, 0.025),
    generations_97.5 = quantile(T_fit * 2 * N_anc, 0.975),
    within_2fold_size = mean(nu_fit / nu_true >= 0.5 & nu_fit / nu_true <= 2),
    within_2fold_time = mean(T_fit / T_true >= 0.5 & T_fit / T_true <= 2),
    .groups = "drop"
  )
write_csv(error_table, file.path(out_dir, "fig3_error_summary.csv"))
print(error_table)


#----------------------------------------------------------
# Shared styling (matches Fig. 1)
#----------------------------------------------------------

theme_fig <- theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey80", color = "grey20"),
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.tag = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank()
  )

log_labels <- scales::trans_format("log10", scales::math_format(10^.x))


#----------------------------------------------------------
# Panel A: population-size histories
# Time on a log scale so the 20-generation event is visible
# even when some replicates infer much older events.
#----------------------------------------------------------

gen_min <- 1
gen_max <- 1e4

history_path <- function(nu, T) {
  t_gen <- pmin(pmax(T * 2 * N_anc, gen_min), gen_max)
  tibble(generations = c(gen_min, t_gen, t_gen, gen_max),
         Ne = c(nu * N_anc, nu * N_anc, N_anc, N_anc))
}

hist_df <- ok %>%
  select(size_long, size_short, replicate, nu_fit, T_fit) %>%
  mutate(path = map2(nu_fit, T_fit, history_path)) %>%
  unnest(path)

median_df <- ok %>%
  group_by(size_long, size_short) %>%
  summarize(nu = median(nu_fit), T = median(T_fit), .groups = "drop") %>%
  mutate(path = map2(nu, T, history_path)) %>%
  unnest(path)

true_df <- history_path(nu_true, T_true)

true_label <- sprintf("True: %s \u2192 %s individuals,\n%s generations ago",
                      scales::comma(N_anc), scales::comma(nu_true * N_anc),
                      scales::comma(T_true * 2 * N_anc))

# Smaller sample drawn first so the larger sample's lines sit on top
replicate_layers <- map(rev(size_short), function(lab) {
  geom_path(data = filter(hist_df, size_short == lab),
            aes(generations, Ne, group = replicate, color = size_short),
            alpha = 0.12, linewidth = 0.5)
})

panel_a <- ggplot() +
  replicate_layers +
  geom_path(data = true_df, aes(generations, Ne, linetype = "True history"),
            color = "black", linewidth = 1.2) +
  geom_path(data = median_df, aes(generations, Ne, color = size_short),
            linewidth = 1.8) +
  annotate("text", x = 1.3, y = 3, label = true_label, hjust = 0, vjust = 0,
           size = 16 / .pt) +
  scale_color_manual(values = size_colors, breaks = size_short, name = "Haplotypes") +
  scale_linetype_manual(values = c("True history" = "dashed"), name = NULL) +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1, linewidth = 1.8)),
    linetype = guide_legend(order = 2, override.aes = list(linewidth = 1.2),
                            keywidth = unit(2.5, "lines"))
  ) +
  scale_x_log10(breaks = 10^(0:4), labels = scales::comma) +
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = scales::comma) +
  coord_cartesian(xlim = c(gen_min, gen_max), ylim = c(1, 2500)) +
  labs(x = "Generations before present (log scale)",
       y = "Effective population size (log scale)") +
  theme_fig


#----------------------------------------------------------
# Panel B: DFE, density per log10 |s| interval
# (gamma on 2Na|s|, transformed to |s|, then to log10 |s|)
#----------------------------------------------------------

s_grid <- 10^seq(-6, 0, length.out = 400)

density_log10 <- function(alpha, beta) {
  dgamma(2 * N_anc * s_grid, shape = alpha, scale = beta) * 2 * N_anc * s_grid * log(10)
}

dfe_df <- ok %>%
  mutate(d = map2(alpha_fit, beta_fit,
                  ~ tibble(s = s_grid, density = density_log10(.x, .y)))) %>%
  select(size_short, d) %>%
  unnest(d) %>%
  group_by(size_short, s) %>%
  summarize(
    median = median(density),
    lower = quantile(density, 0.025),
    upper = quantile(density, 0.975),
    .groups = "drop"
  )

true_dfe <- tibble(s = s_grid, density = density_log10(alpha_true, beta_true))

line_colors <- c("True DFE" = "black", size_colors)

panel_b <- ggplot(dfe_df, aes(s)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = size_short), alpha = 0.18) +
  geom_line(aes(y = median, color = size_short), linewidth = 1.3) +
  geom_line(data = true_dfe, aes(s, density, color = "True DFE"),
            linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(values = line_colors, breaks = "True DFE", name = NULL) +
  scale_fill_manual(values = size_colors, guide = "none") +
  guides(color = guide_legend(override.aes = list(linetype = "dashed"))) +
  scale_x_log10(
    breaks = 10^(-6:0), labels = log_labels,
    sec.axis = sec_axis(~ . * 2 * N_anc, name = expression(2 * N[a] * "|s|"),
                        breaks = 10^(-2:3), labels = log_labels)
  ) +
  labs(x = "Selection coefficient |s| (log scale)",
       y = expression("Density per " * log[10] * "|s| unit")) +
  theme_fig


#----------------------------------------------------------
# Panel C: share of new mutations in selection classes
# Classes use the ancestral Ne (the scale the DFE is inferred on)
#----------------------------------------------------------

class_breaks <- c(0, 1, 10, 100, Inf)
class_labels <- c("Effectively\nneutral\n(< 1)", "Weakly\ndeleterious\n(1-10)",
                  "Moderately\ndeleterious\n(10-100)", "Strongly\ndeleterious\n(> 100)")

class_props <- function(alpha, beta) {
  tibble(class = factor(class_labels, levels = class_labels),
         prop = diff(pgamma(class_breaks, shape = alpha, scale = beta)))
}

class_df <- ok %>%
  mutate(p = map2(alpha_fit, beta_fit, class_props)) %>%
  select(size_short, p) %>%
  unnest(p) %>%
  group_by(size_short, class) %>%
  summarize(
    median = median(prop),
    lower = quantile(prop, 0.025),
    upper = quantile(prop, 0.975),
    .groups = "drop"
  )

true_class <- class_props(alpha_true, beta_true) %>% mutate(label = "True DFE")

# Numbers for the text: true vs inferred share of mutations in each class
class_df %>%
  left_join(select(true_class, class, true = prop), by = "class") %>%
  mutate(class = str_replace_all(class, "\n", " ")) %>%
  select(sample_size = size_short, class, true, median, lower, upper) %>%
  write_csv(file.path(out_dir, "fig3_class_proportions.csv"))

panel_c <- ggplot() +
  geom_col(data = true_class, aes(class, prop, fill = label),
           color = "black", width = 0.75) +
  geom_pointrange(data = class_df,
                  aes(class, median, ymin = lower, ymax = upper, color = size_short),
                  position = position_dodge(width = 0.45), size = 0.7, linewidth = 1.1) +
  scale_fill_manual(values = c("True DFE" = "grey85"), name = NULL) +
  scale_color_manual(values = size_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression("Selection class (" * 2 * N[a] * "|s|)"),
       y = "Share of new mutations") +
  theme_fig


#----------------------------------------------------------
# Assemble and save
#----------------------------------------------------------

# Main figure: A on top, B and C below
fig3 <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A")

fig3
ggsave(file.path(out_dir, "Fig3_demography_DFE.pdf"), fig3,
       width = 16, height = 13, device = cairo_pdf)

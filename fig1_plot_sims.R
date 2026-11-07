library(tidyverse)

#==========================================================
# Fig. 1 plots from the summary CSVs written by fig1_simulate.R
#
# Runs in seconds; no simulation needed. Set csv_dir to the
# folder with the CSVs (the SLURM wrapper sets FIG1_OUTDIR).
# PDFs are saved to the same folder.
#==========================================================

csv_dir <- Sys.getenv("FIG1_OUTDIR", unset = getwd())
setwd(csv_dir)

if (!interactive()) pdf(NULL)  # no stray Rplots.pdf in batch mode

rd <- function(f) read_csv(f, show_col_types = FALSE)


#----------------------------------------------------------
# Labels (rebuilt from the data, so they always match the run)
#----------------------------------------------------------

# Parsed so "ST" is a subscript; quoted so 0.050 prints as "0.050"
fst_label_levels <- function(fsts) paste0('F[ST] == "', sprintf("%.3f", fsts), '"')

add_labels <- function(df) {
  df %>%
    mutate(
      sample_size = factor(n_haplotypes_per_pop,
                           levels = sort(unique(n_haplotypes_per_pop))),
      fst_lab = factor(fst_label_levels(fst),
                       levels = fst_label_levels(sort(unique(fst))))
    )
}

haplotype_labeller <- function(x) paste(x, "haplotypes")


#----------------------------------------------------------
# Shared styling
#----------------------------------------------------------

theme_fig <- theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey80", color = "grey20"),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 16),
    plot.caption = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

size_scales <- list(
  scale_color_viridis_d(end = 0.85, name = "Haplotypes per population"),
  scale_fill_viridis_d(end = 0.85, name = "Haplotypes per population")
)

y_cor_lab <- "|Correlation between PC1 and population|"
annot_size <- 14 / .pt  # 14 pt annotation text
tile_text_size <- 13 / .pt  # 13 pt heatmap cell labels


# Line plot with ribbon and dashed null line (used 3 times)
line_plot <- function(df, y, lower, upper, null_col) {
  ggplot(df, aes(L, {{ y }}, color = sample_size, fill = sample_size, group = sample_size)) +
    geom_ribbon(aes(ymin = {{ lower }}, ymax = {{ upper }}), alpha = 0.12, color = NA) +
    geom_line(aes(y = {{ null_col }}), linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.2) +
    size_scales +
    scale_x_log10(breaks = sort(unique(df$L)), labels = scales::comma) +
    labs(x = "Number of independent SNPs") +
    theme_fig
}

# Detection-power heatmap (used for markers and for LD)
power_heatmap <- function(df, x_var, x_lab, outline_caption) {
  ggplot(df, aes(factor({{ x_var }}), sample_size, fill = power * 100)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_tile(data = filter(df, predicted_detectable),
              fill = NA, color = "black", linewidth = 0.7) +
    geom_text(aes(label = round(power * 100), color = power > 0.6),
              size = tile_text_size) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
    scale_fill_viridis_c(limits = c(0, 100), name = "Detection rate (%)") +
    scale_x_discrete(labels = function(x) scales::comma(as.numeric(x))) +
    labs(x = x_lab, y = "Haplotypes per population", caption = outline_caption) +
    theme_fig +
    theme(panel.grid = element_blank(), legend.key.width = unit(1.5, "cm"))
}


#==========================================================
# Main text
#==========================================================

#----------------------------------------------------------
# Detection-power heatmap (independent markers)
#----------------------------------------------------------

power_df <- rd("fig_detection_power.csv") %>% add_labels()

power_plot <- power_heatmap(
  power_df, L, "Number of independent SNPs",
  expression("Detection: |r| above the 95th percentile of the null. Outlined: " *
               F[ST] %*% sqrt("individuals" %*% "SNPs") > 1)
) +
  facet_wrap(~fst_lab, nrow = 1, labeller = label_parsed)

power_plot
ggsave("Fig_detection_power.pdf", power_plot, width = 20, height = 7)


#----------------------------------------------------------
# |r| vs number of SNPs
#----------------------------------------------------------

cor_df <- rd("fig_correlation.csv") %>% add_labels()

correlation_plot <- line_plot(cor_df, median_correlation, correlation_lower,
                              correlation_upper, null_q95) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~fst_lab, nrow = 1, labeller = label_parsed) +
  labs(
    y = y_cor_lab,
    caption = expression("Dashed lines: 95th percentile of |r| with no structure ("*F[ST] == 0*")")
  )

correlation_plot
ggsave("Fig_PC1_population_correlation.pdf", correlation_plot, width = 20, height = 6.5)


#----------------------------------------------------------
# Collapse plot
#----------------------------------------------------------

collapse_df <- rd("fig_collapse.csv") %>% add_labels()

collapse_plot <- ggplot(collapse_df, aes(bbp_scale, adjusted_correlation, color = sample_size)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "loess", span = 0.5,
              se = FALSE, color = "black", linewidth = 0.7) +
  annotate("text", x = 1, y = 1.04, label = "predicted threshold",
           hjust = -0.05, size = annot_size, color = "grey30") +
  scale_color_viridis_d(end = 0.85, name = "Haplotypes per population") +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  coord_cartesian(ylim = c(-0.1, 1.08)) +
  labs(
    x = expression(F[ST] %*% sqrt("number of individuals" %*% "number of SNPs")),
    y = "PC1-population correlation\n(0 = chance, 1 = perfect)"
  ) +
  theme_fig +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

collapse_plot
ggsave("Fig_collapse.pdf", collapse_plot, width = 10, height = 6.5)


#==========================================================
# Supplement
#==========================================================

#----------------------------------------------------------
# Standardized separation
#----------------------------------------------------------

sep_df <- rd("suppFig_separation.csv") %>% add_labels()

separation_plot <- line_plot(sep_df, median_separation, separation_lower,
                             separation_upper, null_sep_q95) +
  facet_wrap(~fst_lab, nrow = 1, scales = "free_y", labeller = label_parsed) +
  labs(
    y = "Standardized separation along PC1",
    caption = expression("Dashed lines: 95th percentile of separation with no structure ("*F[ST] == 0*")")
  )

separation_plot
ggsave("FigS_standardized_separation.pdf", separation_plot, width = 20, height = 6.5)


#----------------------------------------------------------
# Variance explained by PC1
#----------------------------------------------------------

pve_df <- rd("suppFig_PVE1.csv") %>% add_labels()

pve_plot <- line_plot(pve_df, median_PVE1, PVE1_lower, PVE1_upper, null_PVE1_median) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  facet_wrap(~fst_lab, nrow = 1, labeller = label_parsed) +
  labs(
    y = "Variance explained by PC1",
    caption = expression("Dashed lines: median variance explained with no structure ("*F[ST] == 0*")")
  )

pve_plot
ggsave("FigS_PVE1.pdf", pve_plot, width = 20, height = 6.5)


#----------------------------------------------------------
# LD detection-power heatmap
# Rows = LD design, columns = FST; x = total SNPs genotyped.
# Outlined cells: theory with each linkage block counted as
# one independent locus.
#----------------------------------------------------------

ld_scenario <- function(df) {
  df <- df %>%
    mutate(scenario = ifelse(
      block_size == 1,
      '"Independent SNPs"',
      paste0('"', block_size, ' SNPs/block, " * r^2 == ', r2)
    ))
  levels <- df %>%
    distinct(block_size, r2, scenario) %>%
    arrange(block_size, desc(r2)) %>%
    pull(scenario)
  mutate(df, scenario = factor(scenario, levels = levels))
}

ld_power_df <- rd("suppFig_LD_detection_power.csv") %>% add_labels() %>% ld_scenario()

ld_power_plot <- power_heatmap(
  ld_power_df, L_total, "Total number of SNPs genotyped",
  expression("Outlined: " * F[ST] %*% sqrt("individuals" %*% "linkage blocks") > 1 *
               "  (each block counted as one independent locus)")
) +
  facet_grid(scenario ~ fst_lab, labeller = label_parsed)

ld_power_plot
ggsave("FigS_LD_detection_power.pdf", ld_power_plot, width = 20, height = 24)


#----------------------------------------------------------
# LD |r| line plot
#----------------------------------------------------------

ld_df <- rd("suppFig_LD_correlation.csv") %>%
  add_labels() %>%
  mutate(
    block_lab = factor(
      ifelse(block_size == 1, "1 (independent)", as.character(block_size)),
      levels = c("1 (independent)", as.character(sort(unique(block_size[block_size > 1]))))
    ),
    r2_lab = factor(r2, levels = sort(unique(r2), decreasing = TRUE))
  )

ld_plot <- ggplot(
  ld_df,
  aes(n_blocks, median_correlation, color = block_lab, shape = block_lab,
      linetype = r2_lab, group = interaction(block_lab, r2_lab))
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.3) +
  scale_color_viridis_d(end = 0.85, name = "SNPs per block") +
  scale_shape_manual(values = c(16, 17, 15), name = "SNPs per block") +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted"),
    labels = parse(text = paste0("r^2 == ", levels(ld_df$r2_lab))),
    name = "Within-block LD"
  ) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
  facet_grid(sample_size ~ fst_lab,
             labeller = labeller(fst_lab = label_parsed, sample_size = haplotype_labeller)) +
  labs(x = "Number of linkage blocks (total SNPs / block size)", y = y_cor_lab) +
  theme_fig +
  theme(legend.box = "vertical")

ld_plot
ggsave("FigS_LD_correlation.pdf", ld_plot, width = 20, height = 20)


#----------------------------------------------------------
# Example PC1 vs PC2 plots: one figure per sample size,
# FST down the rows (top = no structure), SNPs across columns
#----------------------------------------------------------

example_df <- rd("suppFig_example_PCA_scores.csv") %>%
  add_labels() %>%
  mutate(L_lab = factor(paste(scales::comma(L), "SNPs"),
                        levels = paste(scales::comma(sort(unique(L))), "SNPs")))

make_example_plot <- function(n_hap) {
  ggplot(filter(example_df, n_haplotypes_per_pop == n_hap), aes(PC1, PC2, color = pop)) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_color_manual(values = c("population 1" = "grey20", "population 2" = "#E69F00"),
                       name = NULL) +
    facet_grid(fst_lab ~ L_lab, scales = "free",
               labeller = labeller(fst_lab = label_parsed)) +
    labs(title = paste(n_hap, "haplotypes per population")) +
    theme_fig +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

for (n_hap in sort(unique(example_df$n_haplotypes_per_pop))) {
  ggsave(sprintf("FigS_example_PCAs_n%d.pdf", n_hap),
         make_example_plot(n_hap), width = 12, height = 18)
}
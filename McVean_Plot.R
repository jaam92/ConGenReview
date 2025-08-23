library(ggplot2)
library(tidyverse)
library(patchwork)
setwd('~/Documents/Perspective_BMC/')
tiger_metadata <- read.csv('per_metadata_tigers.csv')
tiger_metadata2 <- read.csv('per_metadata_tigers2.csv')


# sub10 --------------------------------------------------------------------
pca_vec_100k <- read.table('bengal-malayan-sub10_maf01_thin25kb_100k_fix.eigenvec',head=F)
pca_val_100k <- read.table('bengal-malayan-sub10_maf01_thin25kb_100k_fix.eigenval',head=F)

pca_vec_10k <- read.table('bengal-malayan-sub10_maf01_thin25kb_10k_fix.eigenvec',head=F)
pca_val_10k <- read.table('bengal-malayan-sub10_maf01_thin25kb_10k_fix.eigenval',head=F)

pca_vec_1k <- read.table('bengal-malayan-sub10_maf01_thin25kb_1k_fix.eigenvec',head=F)
pca_val_1k <- read.table('bengal-malayan-sub10_maf01_thin25kb_1k_fix.eigenval',head=F)

pca_vec_100 <- read.table('bengal-malayan-sub10_maf01_thin25kb_100_fix.eigenvec',head=F)
pca_val_100 <- read.table('bengal-malayan-sub10_maf01_thin25kb_100_fix.eigenval',head=F)

pca_vec_10 <- read.table('bengal-malayan-sub10_maf01_thin25kb_10_fix.eigenvec',head=F)
pca_val_10 <- read.table('bengal-malayan-sub10_maf01_thin25kb_10_fix.eigenval',head=F)

pca_val_100k$V2 <- (pca_val_100k$V1/20)*100
pca_val_10k$V2 <- (pca_val_10k$V1/20)*100
pca_val_1k$V2 <- (pca_val_1k$V1/20)*100
pca_val_100$V2 <- (pca_val_100$V1/20)*100
pca_val_10$V2 <- (pca_val_10$V1/20)*100

pca_vec_100k <- pca_vec_100k %>%
  dplyr::rename(Individual = V2)
pca_vec_10k <- pca_vec_10k %>%
  dplyr::rename(Individual = V2)
pca_vec_1k <- pca_vec_1k %>%
  dplyr::rename(Individual = V2)
pca_vec_100 <- pca_vec_100 %>%
  dplyr::rename(Individual = V2)
pca_vec_10 <- pca_vec_10 %>%
  dplyr::rename(Individual = V2)

pca_vec_100klabel <- left_join(pca_vec_100k, tiger_metadata2, by="Individual")
pca_vec_10klabel <- left_join(pca_vec_10k, tiger_metadata2, by="Individual")
pca_vec_1klabel <- left_join(pca_vec_1k, tiger_metadata2, by="Individual")
pca_vec_100label <- left_join(pca_vec_100, tiger_metadata2, by="Individual")
pca_vec_10label <- left_join(pca_vec_10, tiger_metadata2, by="Individual")


pc1_wide <- tibble(
  label = pca_vec_100klabel$Group,   # change if your label column has a different name
  `100k` = pca_vec_100klabel[[3]],
  `10k`  = pca_vec_10klabel[[3]],
  `1k`   = pca_vec_1klabel[[3]],
  `100`  = pca_vec_100label[[3]],
  `10`  = pca_vec_10label[[3]]
)

pc1_long <- pc1_wide %>%
  pivot_longer(cols = c(`100k`, `10k`, `1k`, `100`,`10`),
               names_to = "set",
               values_to = "PC1") %>%
  mutate(n_snps = case_when(
    set == "100k" ~ 100000,
    set == "10k" ~ 10000,
    set == "1k"  ~ 1000,
    set == "100" ~ 100,
    set == "10" ~ 10
  ))


cols <- c("Bengal" = "midnightblue",  # orange
          "Malayan" = "lightskyblue")

a <- ggplot(pc1_long, aes(x = n_snps, y = PC1)) +
  geom_point(aes(color = label), size = 2, alpha = 0.9,
             position = position_jitter(width = 0.04, height = 0)) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000),
                labels = scales::comma) +
  labs(x = "", y = "PC1 (n = 20 haplotypes)") +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",   # put items side-by-side
    legend.background = element_rect(
      fill = scales::alpha("white", 0.9),
      colour = "grey30", linewidth = 0.4
    ),
    legend.key = element_rect(fill = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# sub4 -------------------------------------------------------------------
pca_vec_100k <- read.table('bengal-malayan-sub4_maf01_thin25kb_100k_fix.eigenvec',head=F)
pca_val_100k <- read.table('bengal-malayan-sub4_maf01_thin25kb_100k_fix.eigenval',head=F)

pca_vec_10k <- read.table('bengal-malayan-sub4_maf01_thin25kb_10k_fix.eigenvec',head=F)
pca_val_10k <- read.table('bengal-malayan-sub4_maf01_thin25kb_10k_fix.eigenval',head=F)

pca_vec_1k <- read.table('bengal-malayan-sub4_maf01_thin25kb_1k_fix.eigenvec',head=F)
pca_val_1k <- read.table('bengal-malayan-sub4_maf01_thin25kb_1k_fix.eigenval',head=F)

pca_vec_100 <- read.table('bengal-malayan-sub4_maf01_thin25kb_100_fix.eigenvec',head=F)
pca_val_100 <- read.table('bengal-malayan-sub4_maf01_thin25kb_100_fix.eigenval',head=F)

pca_vec_10 <- read.table('bengal-malayan-sub4_maf01_thin25kb_10_fix.eigenvec',head=F)
pca_val_10 <- read.table('bengal-malayan-sub4_maf01_thin25kb_10_fix.eigenval',head=F)

pca_val_100k$V2 <- (pca_val_100k$V1/8)*100
pca_val_10k$V2 <- (pca_val_10k$V1/8)*100
pca_val_1k$V2 <- (pca_val_1k$V1/8)*100
pca_val_100$V2 <- (pca_val_100$V1/8)*100
pca_val_10$V2 <- (pca_val_10$V1/8)*100

pca_vec_100k <- pca_vec_100k %>%
  dplyr::rename(Individual = V2)
pca_vec_10k <- pca_vec_10k %>%
  dplyr::rename(Individual = V2)
pca_vec_1k <- pca_vec_1k %>%
  dplyr::rename(Individual = V2)
pca_vec_100 <- pca_vec_100 %>%
  dplyr::rename(Individual = V2)
pca_vec_10 <- pca_vec_10 %>%
  dplyr::rename(Individual = V2)

pca_vec_100klabel <- left_join(pca_vec_100k, tiger_metadata2, by="Individual")
pca_vec_10klabel <- left_join(pca_vec_10k, tiger_metadata2, by="Individual")
pca_vec_1klabel <- left_join(pca_vec_1k, tiger_metadata2, by="Individual")
pca_vec_100label <- left_join(pca_vec_100, tiger_metadata2, by="Individual")
pca_vec_10label <- left_join(pca_vec_10, tiger_metadata2, by="Individual")


pc1_wide <- tibble(
  label = pca_vec_100klabel$Group,   # change if your label column has a different name
  `100k` = pca_vec_100klabel[[3]],
  `10k`  = pca_vec_10klabel[[3]],
  `1k`   = pca_vec_1klabel[[3]],
  `100`  = pca_vec_100label[[3]],
  `10`  = pca_vec_10label[[3]]
)

pc1_long <- pc1_wide %>%
  pivot_longer(cols = c(`100k`, `10k`, `1k`, `100`,`10`),
               names_to = "set",
               values_to = "PC1") %>%
  mutate(n_snps = case_when(
    set == "100k" ~ 100000,
    set == "10k" ~ 10000,
    set == "1k"  ~ 1000,
    set == "100" ~ 100,
    set == "10" ~ 10
  ))


cols <- c("Bengal" = "midnightblue",  # orange
          "Malayan" = "lightskyblue")

b <- ggplot(pc1_long, aes(x = n_snps, y = PC1)) +
  geom_point(aes(color = label), size = 2, alpha = 0.9,
             position = position_jitter(width = 0.04, height = 0)) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000),
                labels = scales::comma) +
  labs(x = "Number of independent SNPs", y = "PC1 (n = 8 haplotypes)") +
  theme_bw() +
  theme(legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))


a/b


# add fst if desired ------------------------------------------------------

fst_df <- tribble(
  ~set,   ~n_snps, ~fst,
  "10",      10,   0.081762,
  "100",    100,   0.127732,
  "1k",    1000,   0.11900,
  "10k",  10000,   0.116787,
  "100k",100000,   0.118406
)

fst_df <- fst_df %>% semi_join(pc1_long %>% distinct(n_snps), by = "n_snps")

# y-position a bit above the tallest point in each bin
tops <- pc1_long %>% group_by(n_snps) %>% summarise(top = max(PC1, na.rm = TRUE))
ann  <- fst_df %>%
  left_join(tops, by = "n_snps") %>%
  mutate(y = top + 0.035,                              # vertical offset
         label = paste0("italic(F)[ST]==", sprintf('%.3f', fst)))  # FST with F italic & ST subscript

y_needed <- max(ann$y, na.rm = TRUE) + 0.02            # make room for labels

label_map <- c("Bengal" = "Bengal (n = 25)",
               "Malayan" = "Malayan (n = 24)")

a <- ggplot(pc1_long, aes(x = n_snps, y = PC1)) +
  geom_point(aes(color = label), size = 2, alpha = 0.9,
             position = position_jitter(width = 0.04, height = 0)) +
  scale_color_manual(
    values = cols, name = NULL,
    labels = function(b) { out <- b
    out[out == "Bengal"]  <- "Bengal (n = 25)"
    out[out == "Malayan"] <- "Malayan (n = 24)"
    out
    }
  ) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000),
                labels = scales::comma) +
  labs(x = "Number of independent SNPs", y = "PC1") +
  theme_bw() +
  theme(
    legend.position = c(0.75, 0.15),      # inside plot: (x,y) in [0,1]
    legend.justification = c(0, 1),       # anchor the top-left corner
    legend.background = element_rect(
      fill = scales::alpha("white", 0.9), # boxed legend
      colour = "grey30", linewidth = 0.4
    ),
    legend.key = element_rect(fill = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
geom_label(data = ann,
           aes(x = n_snps, y = y, label = label),
           parse = TRUE, fill = "white", color = "grey30",
           label.size = 0.3, label.padding = unit(0.12, "lines"),
           size = 3.3) +
  expand_limits(y = y_needed)





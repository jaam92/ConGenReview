#Load library
library(tidyverse)

#function to generate genotypes from each population where we specificy number of individuals (n_per_pop), markers (L), and FST
simulate_data <- function(n_per_pop = 5, L = 1000, fst = 0.01) {
  p0 <- runif(L, 0.1, 0.9)
  delta <- sqrt(fst * p0 * (1 - p0)) #make sure FST is approx. the intended value
  p1 <- p0 + delta
  p2 <- p0 - delta
  p1[p1 > 1] <- 1 #adjust if the delta coefficient made allele freq. go over 1
  p2[p2 < 0] <- 0 #adjust if the delta coefficient made allele freq. less than 0
  
  pop1 <- replicate(n_per_pop, rbinom(L, 2, p1)) #this gives diploid GTs
  pop2 <- replicate(n_per_pop, rbinom(L, 2, p2)) #this gives diploid GTs
  
  dat <- data.frame(t(cbind(pop1, pop2)))
  dat$pop <- rep(c("population 1", "population 2"), each = n_per_pop)
  dat
}

#vary L and FST and number of diploid samples
Ls  <- c(10, 100, 1000, 10000, 100000)
FSTs <- seq(0.01, 0.10, by = 0.01)
n_vals <- c(10, 4) 
set.seed(120)
results <- list()

#Loop through everything
for (fst in FSTs) {
  for (L in Ls) {
    for (n_per_pop in n_vals) {
      
      dat <- simulate_data(n_per_pop = n_per_pop, L = L, fst = fst)
      X <- dat %>% 
        select(-pop)
      
      keep <- apply(X, 2, var) > 0
      X <- X[, keep]
      
      pca <- prcomp(X, scale. = TRUE)
      pve1 <- pca$sdev[1]^2 / sum(pca$sdev^2) #proportion of variance explained by pc1
      
      df <- data.frame(
        PC1 = pca$x[,1],
        PVE1 = pve1,
        pop = dat$pop,
        L = L,
        fst = fst,
        n_per_pop = n_per_pop
      )
      
      results[[paste(L, fst, n_per_pop, sep = "_")]] <- df
    }
  }
}

#Put everything together
df_all <- bind_rows(results)


# Make facet label variables
df_all <- df_all %>%
  mutate(
    fst_lab = paste0("F[ST]==", fst),
    n_lab   = paste0("n==", 2*n_per_pop) #number of haplotypes is 2*num dips
  )

#Choose colors
cols <- c(
  "population 1" = "#67271BFF",
  "population 2" = "#E7A79BFF"
)

#Plot PC1
ggplot(df_all, aes(x = L, y = PC1)) +
  geom_point(
    aes(color = pop),
    size = 2,
    alpha = 0.9,
    position = position_jitter(width = 0.04, height = 0)
  ) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000, 100000),
    labels = scales::comma
  ) +
  labs(
    x = "Number of independent SNPs",
    y = "PC1"
  ) +
  facet_grid(
    n_lab ~ fst_lab,   # rows: sample size, cols: FST
    labeller = label_parsed
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
    legend.text = element_text(size = 18),
  )

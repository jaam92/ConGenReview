library(tidyverse)
library(magick)
library(pdftools)
library(grid)

library(patchwork)

setwd("C:/Users/tferrari/Desktop/Professional/My Papers/BMC_congen_review")

fig3_image <- image_read_pdf("fig3_error.pdf", pages = 1)
fig3 <- wrap_elements(rasterGrob(fig3_image))

fig4_image <- image_read_pdf("fig4_diffs_barplot.pdf", pages = 1)
fig4 <- wrap_elements(rasterGrob(fig4_image)) +
  labs(title = "C.") +
  theme(plot.title = element_text(hjust = 0.2, face = "bold", size = 16))

new_fig3 <- fig3 / fig4 + plot_layout(heights = c(2, 2.5))
ggsave(new_fig3, 
       filename = "new_fig3_dfe_error.pdf",
       device = cairo_pdf, width = 10, height = 10, units = "in")

library(tidyverse)

fig_size <- 3

melt <- function(mat) {
  mat %>%
    magrittr::set_colnames(NULL) %>%
    magrittr::set_rownames(NULL) %>%
    reshape2::melt(varnames = c("y", "x"), value.name = "value") %>%
    as_tibble()
}
raster <- function(meltdf, legend = FALSE, palette = "RdYlBu") {
  g <-
    ggplot(meltdf) +
    geom_raster(aes(x, y, fill = value)) +
    scale_fill_distiller(palette = palette, breaks = c(0, 1), labels = c("Min", "Max")) +
    dynplot::theme_graph() +
    labs(fill = "Expression") +
    coord_cartesian(expand = FALSE)
  if (!legend) {
    g <- g + theme(legend.position = "none")
  }
  g
}


dat <- dyntoy::generate_dataset(id = "defence", model = "multifurcating", num_cells = 100, num_features = 101)

counts <-
  dat$counts %>%
  as.matrix()
colnames(counts) <- gsub("G", "gene ", colnames(counts))
rownames(counts) <- gsub("C", "cell ", rownames(counts))
write_tsv(data.frame(counts, check.names = FALSE) %>% rownames_to_column("cell name"), "figures/4_singlecellomics/rnaseq_matrix.tsv")

g <-
  counts %>%
  melt() %>%
  raster()

ggsave("figures/4_singlecellomics/rnaseq_plot.pdf", g, width = 3, height = 3)


countsmelt <-
  reshape2::melt(counts, varnames = c("cell", "gene"), value.name = "count")
countsmelt_sel <-
  countsmelt %>%
  filter(cell == "cell 1", gene %in% paste0("gene ", 1:5))
g <- ggplot(countsmelt_sel) + geom_bar(aes(gene, count), stat = "identity") + theme_classic() + labs(x = NULL, y = "Aantal RNA transcripten") + coord_flip()
ggsave("figures/4_singlecellomics/rnaseq_2dplot.pdf", g, width = 3, height = 3)


g <-
  counts %>%
  melt() %>%
  filter(y == 1, x %in% 1:5) %>%
  raster()

ggsave("figures/4_singlecellomics/rnaseq_heatmap1cell.pdf", g, width = 3, height = 3/5)

scale_fill_distiller(palette = palette, breaks = c(0, 1), labels = c("Min", "Max")) +

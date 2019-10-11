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


# Generate toy datasets ---------------------------------------------------
set.seed(1)

# linear
linear <- dyntoy::generate_dataset(
  model = "linear",
  num_cells = 60,
  num_features = 40,
  dropout_probability_factor = 0
)

expr_li <- as.matrix(linear$expression) %>% dynutils::scale_quantile() %>% t() %>% dynutils::scale_quantile()

dropout_weight <- 1 / (as.vector(expr_li)+1)
rm_ix <- sample(seq_along(dropout_weight), length(dropout_weight) * .2, prob = dropout_weight, replace = FALSE)

expr_li_nodrop <- expr_li
expr_li[rm_ix] <- 0

ordr_li <- hclust(dist(expr_li))$order
ordc_li <- hclust(dist(t(expr_li)))$order

# branching
set.seed(1)
branching <- dyntoy::generate_dataset(
  model = dyntoy::model_binary_tree(num_branchpoints = 2L),
  num_cells = 100,
  num_features = 100,
  allow_tented_progressions = FALSE,
  sample_mean_count = function() rnorm(1, 1000, 100),
  sample_dispersion_count = function(mean) map_dbl(mean, ~runif(1, ./10, ./4)),
  dropout_probability_factor = 10
)
expr_br <- as.matrix(branching$expression)
expr_br_sc <- expr_br %>% dynutils::scale_quantile() %>% t() %>% dynutils::scale_quantile()

ordr_br <- hclust(dist(expr_br_sc))$order
ordc_br <- hclust(dist(t(expr_br_sc)))$order

# Experiment --------------------------------------------------------------
ggsave(
  "fig/comptools/1_experiment_expression.pdf",
  expr_li %>% melt() %>% raster(),
  width = fig_size,
  height = fig_size
)
ggsave(
  "fig/comptools/1_legend.pdf",
  expr_li %>% melt() %>% raster(legend = TRUE),
  width = fig_size,
  height = fig_size
)

# Imputation --------------------------------------------------------------
ggsave(
  "fig/comptools/2_imputation_input.pdf",
  expr_li[ordr_li, ordc_li] %>% melt() %>% raster(),
  width = fig_size,
  height = fig_size
)

ggsave(
  "fig/comptools/2_imputation_output.pdf",
  expr_li_nodrop[ordr_li, ordc_li] %>% melt() %>% raster(),
  width = fig_size,
  height = fig_size
)

# Integration -------------------------------------------------------------
g1 <- expr_li[sample(ordr_li), ordc_li[1:20] %>% sample()] %>% melt() %>% raster(palette = "RdYlGn")
g2 <- expr_li[sample(ordr_li), ordc_li[21:40] %>% sample()] %>% melt() %>% raster(palette = "RdBu")
g3 <- expr_li[sample(ordr_li), ordc_li[41:60] %>% sample()] %>% melt() %>% raster(palette = "PiYG")
g <- patchwork::wrap_plots(g1, g2, g3, nrow = 1)

ggsave(
  "fig/comptools/3_integration_input.pdf",
  g,
  width = fig_size,
  height = fig_size
)

g1 <- expr_li[ordr_li, ordc_li[1:20]] %>% melt() %>% raster(palette = "RdYlGn")
g2 <- expr_li[ordr_li, ordc_li[21:40]] %>% melt() %>% raster(palette = "RdBu")
g3 <- expr_li[ordr_li, ordc_li[41:60]] %>% melt() %>% raster(palette = "PiYG")
g <- patchwork::wrap_plots(g1, g2, g3, nrow = 1)

ggsave(
  "fig/comptools/3_integration_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)


# Dimensionality reduction ------------------------------------------------
ggsave(
  "fig/comptools/4_dimred_input.pdf",
  expr_br_sc[ordr_br, ordc_br] %>% melt() %>% raster(),
  width = fig_size,
  height = fig_size
)

dimred_br <- dyndimred::dimred_tsne(expr_br, perplexity = 4)

g <-
  dynplot::plot_dimred(
    branching,
    dimred = dimred_br,
    color_cells = "none",
    plot_trajectory = FALSE,
    size_cells = 1
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()

ggsave(
  "fig/comptools/4_dimred_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)

# Clustering --------------------------------------------------------------



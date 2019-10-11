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
set.seed(2)
branching <- dyntoy::generate_dataset(
  model = dyntoy::model_binary_tree(num_branchpoints = 1L),
  num_cells = 100,
  num_features = 100,
  allow_tented_progressions = FALSE,
  sample_mean_count = function() rnorm(1, 1000, 100),
  sample_dispersion_count = function(mean) map_dbl(mean, ~runif(1, ./10, ./4)),
  dropout_probability_factor = 10
)
expr_br <- as.matrix(branching$expression)
expr_br_sc <- expr_br %>% dynutils::scale_quantile() %>% t() %>% dynutils::scale_quantile()

hclr_br <- hclust(dist(expr_br_sc))
ordr_br <- hclr_br$order
hclc_br <- hclust(dist(t(expr_br_sc)))
ordc_br <- hclc_br$order

# Experiment --------------------------------------------------------------
set.seed(1)
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
    size_cells = 2
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
ggsave(
  "fig/comptools/5_cluster_input.pdf",
  expr_br_sc[ordr_br, ordc_br] %>% melt() %>% raster(),
  width = fig_size,
  height = fig_size
)

cluscut_br <- cutree(hclc_br, k = 3)

clus_ord <- unique(cluscut_br[ordc_br])
rasters <- map(clus_ord, function(cli) {
  ix <- cluscut_br[ordc_br] == cli
  expr_br_sc[ordr_br, ordc_br[ix]] %>% melt() %>% raster()
})
widths <- map_int(clus_ord, ~ sum(. == cluscut_br))
g <- patchwork::wrap_plots(rasters, nrow = 1, widths = widths)

ggsave(
  "fig/comptools/5_cluster_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)

g <-
  dynplot::plot_dimred(
    branching,
    dimred = dimred_br,
    color_cells = "grouping",
    plot_trajectory = FALSE,
    grouping = as.character(cluscut_br),
    size_cells = 2
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()

ggsave(
  "fig/comptools/5_clusterdr_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)


# Trajectory inference ----------------------------------------------------

df <-
  tibble(
    edge = dynwrap::group_onto_trajectory_edges(branching),
    progression = branching$progressions %>% slice(match(colnames(expr_br_sc), cell_id)) %>% pull(percentage),
    orig_index = seq_along(edge)
  ) %>%
  arrange(edge, progression)


edge_ord <- unique(df$edge)
rasters <- map(edge_ord, function(ed) {
  ix <- df %>% filter(edge == ed) %>% pull(orig_index)
  expr_br_sc[ordr_br, ix] %>% melt() %>% raster()
})
widths <- df %>% group_by(edge) %>% summarise(n = n()) %>% pull(n)
g <- patchwork::wrap_plots(rasters, nrow = 1, widths = widths)

ggsave(
  "fig/comptools/6_ti_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)


dimred_br2 <- dyndimred::dimred_mds(expr_br)

g <-
  dynplot::plot_dimred(
    branching,
    dimred = dimred_br2,
    color_cells = "none",
    plot_trajectory = FALSE,
    size_cells = 2
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()

ggsave(
  "fig/comptools/6_tidr_input.pdf",
  g,
  width = fig_size,
  height = fig_size
)

g <-
  dynplot::plot_dimred(
    branching,
    dimred = dimred_br2,
    # color_cells = "none",
    plot_trajectory = TRUE,
    size_cells = 2
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()

ggsave(
  "fig/comptools/6_tidr_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)


# Trajectory alignment ----------------------------------------------------
set.seed(1)


library(dyngen)
backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "doublerep1", num_modules = 6),
  bblego_linear("B", "C", type = "doublerep2", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)
model <-
  initialise_model(
    num_tfs = nrow(backbone$module_info) * 2,
    num_targets = 100 - (nrow(backbone$module_info) * 2),
    num_hks = 0,
    num_cells = 120,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(num_simulations = 10, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network()


model1 <- model
model1$feature_info <- model1$feature_info %>%
  mutate_at(c("ba", "ind"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  mutate_at(c("ba", "ind"), ~ pmin(., 1))
model1$feature_network <- model1$feature_network %>%
  mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
model1 <- model1 %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()

model2 <- model
model2$feature_info <- model2$feature_info %>%
  mutate_at(c("ba", "ind"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  mutate_at(c("ba", "ind"), ~ pmin(., 1))
model2$feature_network <- model2$feature_network %>%
  mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
model2 <- model2 %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()


model12 <- model
model12$feature_info <- model$feature_info %>%
  mutate(
    w = paste0("w_", feature_id),
    x = paste0("x_", feature_id),
    y = paste0("y_", feature_id)
  )
m1gs <- model1$gold_standard
m2gs <- model2$gold_standard
model12$gold_standard <- list(
  meta = bind_rows(
    m1gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("left_", .)),
    m2gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("right_", .))
  ),
  counts = rbind(
    m1gs$counts,
    m2gs$counts
  ),
  network = bind_rows(
    m1gs$network %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
    m2gs$network %>% mutate_at(c("from", "to"), ~paste0("right_", .))
  )
)
m1sim <- model1$simulations
m2sim <- model2$simulations
num_m1sim <- max(m1sim$meta$simulation_i)
model12$simulations <- list(
  meta = bind_rows(
    m1sim$meta %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
    m2sim$meta %>% mutate_at(c("from", "to"), ~paste0("right_", .)) %>% mutate(simulation_i = simulation_i + num_m1sim)
  ),
  counts = rbind(
    m1sim$counts,
    m2sim$counts
  ),
  regulation = rbind(
    m1sim$regulation,
    m2sim$regulation
  ),
  reaction_firings = rbind(
    m1sim$reaction_firings,
    m2sim$reaction_firings
  ),
  reaction_propensities = rbind(
    m1sim$reaction_propensities,
    m2sim$reaction_propensities
  )
)
model12 <- model12 %>%
  dyngen:::calculate_dimred(dimred_premrna = FALSE) %>%
  generate_experiment()

linear12 <- wrap_dataset(model12) %>%
  dynwrap::simplify_trajectory()

dimred12 <- dyndimred::dimred_mds(linear12$expression)

lpt <- linear12$progressions %>% filter(grepl("left_", from)) %>% arrange(percentage) %>% pull(cell_id)
rpt <- linear12$progressions %>% filter(grepl("right_", from)) %>% arrange(percentage) %>% pull(cell_id)

alignment <- dtw::dtw(
  linear12$expression[lpt, ], linear12$expression[rpt, ],
  dist.method = "correlation",
  keep = TRUE
)
al1 <- lpt[alignment$index1s]
al2 <- rpt[alignment$index2s]


df <- data.frame(
  from = al1,
  to = al2,
  from = dimred12[al1,],
  to = dimred12[al2,]
) %>% as_tibble()


g <-
  dynplot::plot_dimred(
    linear12,
    dimred = dimred12,
    plot_trajectory = TRUE,
    size_cells = 2
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()


ggsave(
  "fig/comptools/7_align_input.pdf",
  g,
  width = fig_size,
  height = fig_size
)


g <-
  dynplot::plot_dimred(
    linear12,
    dimred = dimred12,
    plot_trajectory = TRUE,
    size_cells = 2
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian() +
  geom_segment(aes(x = from.comp_1, xend = to.comp_1, y = from.comp_2, yend = to.comp_2), df, colour = "lightgray")


ggsave(
  "fig/comptools/7_align_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)



# Network inference -------------------------------------------------------
set.seed(1)

library(dyngen)
backbone <- backbone_bifurcating()
model_bifur <-
  initialise_model(
    num_tfs = nrow(backbone$module_info) * 2,
    num_targets = 100 - (nrow(backbone$module_info) * 2),
    num_hks = 0,
    num_cells = 120,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(num_simulations = 10, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

bifur <- wrap_dataset(model_bifur, store_dimred = TRUE)

# dynplot::plot_dimred(bifur)

expr_bifur <- bifur$expression %>% as.matrix
expr_bifur_sc <- expr_bifur %>% dynutils::scale_quantile() %>% t() %>% dynutils::scale_quantile()

ordr_li <- hclust(dist(expr_bifur_sc))$order
ordc_li <- hclust(dist(t(expr_bifur_sc)))$order

ggsave(
  "fig/comptools/8_ni_input.pdf",
  expr_bifur_sc[ordr_li, ordc_li] %>% melt() %>% raster(),
  width = fig_size,
  height = fig_size
)

expr_bifur <- expr_bifur[, colMeans(expr_bifur) != 0]
sigmoid_mean <- mean(expr_bifur[expr_bifur != 0])
sigmoid_sd <- sd(expr_bifur[expr_bifur != 0])
targets <- colnames(expr_bifur)
regulators <- model_bifur$feature_info %>% filter(is_tf) %>% pull(feature_id) %>% intersect(targets)
samples <- rownames(expr_bifur)

imps <- pbapply::pblapply(
  seq_along(targets),
  function(i) {
    # cat(i, "\n", sep = "")
    target <- targets[[i]]
    regs <- setdiff(regulators, target)
    target_expr <- scale(expr_bifur[,target])[,1]

    data <- data.frame(
      PREDICT = target_expr,
      as.matrix(expr_bifur[,regs]),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    rf2 <- rangercase::ranger(
      data = data,
      dependent.variable.name = "PREDICT",
      num.threads = 1,
      importance = "permutation",
      num.trees = 1000,
      min.node.size = 5,
      local.importance = TRUE
    )

    imp <- tibble(
      feature_id = names(rf2$variable.importance),
      importance = rf2$variable.importance,
      effect = ifelse(importance == 0, 0, rf2$variable.importance.cor),
      ix = seq_along(feature_id)
    ) %>%
      arrange(desc(importance))

    limp <- Matrix::t(rf2$variable.importance.casewise[imp$feature_id, , drop = FALSE])
    expr_reg_sc <- stats::pnorm(as.matrix(data[rownames(limp),colnames(limp)]), mean = sigmoid_mean, sd = sigmoid_sd)
    limp_sc <- limp * expr_reg_sc

    limp_df <-
      left_join(
        limp %>%
          reshape2::melt(varnames = c("cell_id", "regulator"), value.name = "importance"),
        limp_sc %>%
          reshape2::melt(varnames = c("cell_id", "regulator"), value.name = "importance_sc"),
        by = c("cell_id", "regulator")
      ) %>%
      as_tibble()

    lst(
      importance = imp %>% transmute(
        regulator = factor(feature_id, levels = regulators),
        target = factor(target, levels = targets),
        importance,
        effect
      ),
      importance_sc = limp_df %>% transmute(
        cell_id = factor(as.character(cell_id), levels = samples),
        regulator = factor(as.character(regulator), levels = regulators),
        target = factor(target, levels = targets),
        importance,
        importance_sc
      )
    )
  }
)

importance <- imps %>%
  map_df("importance") %>%
  arrange(desc(importance)) %>%
  head(500) %>%
  mutate_at(c("regulator", "target"), as.character)

feature_info <-
  tibble(
    feature_id = targets,
    is_tf = targets %in% regulators,
    color_by = ifelse(is_tf, "Regulator", "Target")
  ) %>%
  filter(feature_id %in% c(importance$regulator, importance$target))

color_legend <- c(
  "Regulator" = "black",
  "Target" = "darkgray"
)


library(tidygraph)
library(ggraph)
gr <- tbl_graph(nodes = feature_info, edges = importance) %>%
  tidygraph::to_minimum_spanning_tree() %>%
  .$mst

cap <- circle(1, "mm")
arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(2, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(2, "mm"))

g <- ggraph(gr, layout = "kk") +
  geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0 & from != to, colour = "Promotes"), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect < 0, colour = "Represses"), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_node_point(size = 1.5) +
  scale_edge_colour_manual(values = c("Promotes" = "#31a354", "Represses" = "#de2d26")) +
  ggraph::theme_graph(base_family = "Helvetica") +
  coord_equal() +
  theme(legend.position = "bottom") +
  labs(edge_colour = "Interaction type")

ggsave(
  "fig/comptools/8_ni_output.pdf",
  g,
  width = fig_size*1.5,
  height = fig_size*1.5
)


# RNA velocity ------------------------------------------------------------
set.seed(1)

dataset <- dyntoy::generate_dataset(model = "bifurcating", num_cells = 1000, add_prior_information = F, add_velocity = T, allow_tented_progressions = FALSE)
dataset <- dataset %>% add_dimred(dimred = dyndimred::dimred_landmark_mds) %>% add_root()

library(dynplot2)
g <-
  dynplot(dataset) +
  geom_cell_point(color = "grey80", size = 1) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()

ggsave(
  "fig/comptools/9_rnavelocity_input.pdf",
  g,
  width = fig_size,
  height = fig_size
)

g <-
  dynplot(dataset) +
  geom_cell_point(color = "grey80", size = 1) +
  geom_velocity_arrow(stat = stat_velocity_grid(grid_n = 20), arrow = grid::arrow(length = grid::unit(1, "mm"))) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2") +
  coord_cartesian()

ggsave(
  "fig/comptools/9_rnavelocity_output.pdf",
  g,
  width = fig_size,
  height = fig_size
)


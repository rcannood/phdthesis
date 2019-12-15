library(tidyverse)
library(dyngen)

set.seed(3)

backbone <- backbone_bifurcating()

model <-
  initialise_model(
    num_cells = 300,
    num_tfs = 50,
    num_targets = 200,
    num_hks = 0,
    distance_metric = "pearson",
    backbone = backbone,
    tf_network_params = tf_network_default(min_tfs_per_module = 1, sample_num_regulators = function() 2),
    feature_network_params = feature_network_default(target_resampling = 5000),
    kinetics_params = kinetics_default(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(
      census_interval = .1,
      burn_time = 3,
      total_time = 10,
      experiment_params = bind_rows(
        simulation_type_wild_type(num_simulations = 40),
        simulation_type_knockdown(num_simulations = 0)
      )
    ),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8
  )

model <- model %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()

model <- model %>%
  generate_gold_standard()

plot_gold_simulations(model)

model <- model %>% generate_cells()

plot_simulation_expression(model, 1, what = "x")
plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")

simcounts <- model$simulations$counts
dr <- dyndimred::dimred_landmark_mds(simcounts, distance_method = "spearman")
dim(dr)

plot(dr)

drdf <-
  bind_cols(
    model$simulations$meta,
    as.data.frame(dr)
  )

g <-
  ggplot(drdf %>% filter(sim_time >= 0)) +
  geom_path(aes(comp_1, comp_2, group = simulation_i), function(df) df %>% select(-sim_time), colour = "darkgray") +
  geom_point(aes(comp_1, comp_2), colour = "red") +
  gganimate::transition_time(sim_time) +
  theme_classic() +
  labs(x = "Component 1", y = "Component 2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)
gganimate::anim_save("figures/7_dyngen/insilico.gif", g, width = 800, height = 800)

file.remove("figures/7_dyngen/insilico.avi"); system("ffmpeg -f gif -i figures/7_dyngen/insilico.gif -crf 1 -b:v 1M figures/7_dyngen/insilico.avi")
system("convert 'figures/7_dyngen/insilico.gif[0]' figures/7_dyngen/insilico_init.png")

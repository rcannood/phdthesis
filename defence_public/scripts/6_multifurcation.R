library(tidyverse)
library(dyno)

set.seed(6)
dat <- dyntoy::generate_dataset(id = "defence", model = dyntoy::model_multifurcating(num_branchpoints = 2, max_degree = 3), num_cells = 200, num_features = 101, allow_tented_progressions = FALSE, add_velocity = FALSE)

dimred <- dyndimred::dimred_mds(dat$expression) %>% dynutils::scale_uniform()
plot_dimred(dat, dimred = dimred)

gro <- setNames(rep("NA", length(dat$cell_ids)), dat$cell_ids)
g1 <-
  dynplot::plot_dimred(dat, dimred = dimred, plot_trajectory = FALSE, color_cells = "grouping", grouping = gro, size_cells = 1.5, size_milestones = 4) +
  theme_classic() + labs(x = "Component 1", y = "Component 2") +
  scale_x_continuous(breaks = NULL, limits = c(-.5, .5)) + scale_y_continuous(breaks = NULL, limits = c(-.5, .5)) +
  scale_colour_manual(values = c("NA" = "darkgray")) + theme(legend.position = "none")

methods <- list(
  "Monocle DDRTree" = ti_monocle_ddrtree(),
  "PAGA" = ti_paga_tree(),
  "Slingshot" = ti_slingshot(),
  "CellTrails" = ti_celltrails(),
  "cellTree" = ti_celltree_maptpx(),
  "DPT" = ti_dpt()
)

traj_outs <- map(names(methods), function(method_name) {
  method <- methods[[method_name]]
  traj <- infer_trajectory(dat, method)
  plot <- dynplot::plot_dimred(traj, dimred = dimred, color_cells = "grouping", grouping = gro, size_cells = 1.5, size_milestones = 4) +
    labs(title = method_name) +
    scale_colour_manual(values = c("NA" = "lightgray")) + theme(legend.position = "none")
  lst(traj, plot)
})
gs <- patchwork::wrap_plots(map(traj_outs, "plot"), nrow = 2)


g2 <-
  dynplot::plot_dimred(dat, dimred = dimred, color_cells = "grouping", grouping = gro, size_cells = 1.5, size_milestones = 4) + theme_classic() + labs(x = "Component 1", y = "Component 2") +
  scale_x_continuous(breaks = NULL, limits = c(-.5, .5)) + scale_y_continuous(breaks = NULL, limits = c(-.5, .5)) +
  scale_colour_manual(values = c("NA" = "darkgray")) + theme(legend.position = "none")

ggsave("figures/6_multifurcation/plot_unlabelled.pdf", g1, width = 3, height = 3)
ggsave("figures/6_multifurcation/plot_predictions.pdf", gs, width = 3*3, height = 2*2.5)
ggsave("figures/6_multifurcation/plot_labelled.pdf", g2, width = 3, height = 3)

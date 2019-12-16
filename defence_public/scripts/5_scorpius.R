library(SCORPIUS)
library(tidyverse)
library(gganimate)

set.seed(1)

data <- SCORPIUS::ginhoux
expr <- data$expression
sample_info <- data$sample_info

dimred <- SCORPIUS::reduce_dimensionality(expr, ndim = 2)

g2 <-
  draw_trajectory_plot(dimred, point_size = 1) +
  scale_x_continuous(breaks = -Inf) +
  scale_y_continuous(breaks = -Inf)
ggsave("figures/5_scorpius/dimred.pdf", g2, width = 3, height = 3)

its <- 0:10
walk(its, function(it) {
  set.seed(1)
  traj <- SCORPIUS::infer_trajectory(dimred, maxit = it, k = 0)

  pl <- g2 +
    geom_path(aes(Comp1, Comp2), as.data.frame(traj$path), colour = "red") +
    labs(title = paste0("Iteration ", it))
  # ggsave(paste0("figures/5_scorpius/princurve_it", it, ".pdf"), pl, width = 3, height = 3)
  ggsave(paste0("figures/5_scorpius/princurve_it", it, ".png"), pl, width = 3, height = 3)
})

system("convert figures/5_scorpius/princurve_it*.png -delay 1000 figures/5_scorpius/princurve_loop.gif")
file.remove("figures/5_scorpius/princurve_loop.avi"); system("ffmpeg -f image2 -r 1 -i figures/5_scorpius/princurve_it%d.png figures/5_scorpius/princurve_loop.avi")

set.seed(1)
traj <- SCORPIUS::infer_trajectory(dimred, maxit = 10, k = 0)

ix <- sample.int(nrow(expr), 50)
gimp <- SCORPIUS::gene_importances(expr, traj$time)

expr_df <- data.frame(cell = names(traj$time),time = traj$time, as.matrix(expr[,c("Cd74", "Mpo", "Prtn3", "Siglech")])) %>%
  gather(gene, value, -time, -cell)
g6 <-
  ggplot(expr_df) +
  geom_point(aes(time, value), size = 1) +
  geom_smooth(aes(time, value), colour = "red", fill = NA, size = 1) +
  facet_wrap(~gene, scales = "free") +
  theme_classic() +
  scale_x_continuous(breaks = -Inf) +
  scale_y_continuous(breaks = -Inf) +
  theme(strip.background = element_blank(), legend.position = "none") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Pseudotime", y = "Expression")
ggsave("figures/5_scorpius/expr_time.pdf", g6, width = 3, height = 3)

expr_sc <- scale_quantile(expr[,gimp$gene[1:50]])
modules <- extract_modules(expr_sc, traj$time)
draw_trajectory_heatmap(
  expr_sc[ix,], traj$time[ix],
  modules = modules, scale_features = FALSE, border_color = NA,
  annotation_legend = FALSE,
  filename = "figures/5_scorpius//heatmap.pdf",
  width = 4,
  height = 3
)

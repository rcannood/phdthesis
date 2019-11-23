library(tidyverse)
library(scales)

set.seed(1)
df <- tibble(
  x = round(10^runif(6, .5, )),
  n = paste0("Molecule ", seq_along(x)) %>% forcats::fct_rev()
)

ggplot(df) +
  geom_bar(aes(n, x), stat = "identity") +
  theme_classic() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = 10^c(0, 5)) +
  labs(x = NULL, y = "Abundance") +
  coord_flip()
ggsave("defence_internal/figures/1_part_abundance.pdf", width = 4, height = 3)



df <- tribble(
  ~x, ~y, ~name, ~hjust, ~vjust,
  0, 0, "§2: dyngen", .5, 1,
  .4, 1, "§3: dynbenchmark", 0, .5,
  .6, 1, "§8: self-assessment\nin trajectory inference", 0, .5,
  .8, 1, "§5: SCORPIUS", 0, .5,
  1, 1, "§4: dyno", 0, .5,
  1.8, 1, "§6: bred", 1, .5,
  2, 1, "§7: incgraph", 1, .5,
  3, 1, "§9: guidelines\nfor benchmarking", .5, 0
)
titdf <- tibble(
  x = 0:2,
  y = c(1.6, 1.8, 1.8),
  # y = 1.6,
  # y = c(1.6, 1.5, 1.5),
  name = c("Benchmarking", "Trajectory\ninference", "Network\ninference"),
  short = c("Benchmarking", "TI", "NI"),
  col = RColorBrewer::brewer.pal(3, "Set1"),
  hjust = c(.5, 0, 1)
)

rgbs <- grDevices::col2rgb(titdf$col)
mix_cols <- function(df) {
  rgb <-
    map2(df$x, df$y, function(x, y) {
      case_when(
        y == 0 ~ c(1, 1, 1) / 3,
        x <= 1 ~ c(1 - x, x, 0),
        x <= 2 ~ c(0, 2 - x, x - 1),
        TRUE ~ c(x - 2, 0, 3 - x)
      )
    }) %>% do.call(rbind, .)
  rgb2 <- round(rgbs %*% t(rgb))
  df$col <- grDevices::rgb(rgb2[1,], rgb2[2,], rgb2[3,], maxColorValue = 255)
  df
}

df <- df %>% mix_cols()

arrdf <- tribble(
  ~gr, ~x, ~y,
  2, .4, 1,
  2, .76, 1,
  3, 1, 1,
  3, 1.575, 1,
  4, 2, 1,
  4, 2.5, 1
)

g <- ggplot(df) +
  geom_segment(aes(x = x, xend = x, y = 0, yend = 1.05), titdf, colour = "lightgray", linetype = "dashed") +
  geom_path(aes(x, y, group = gr), arrdf, arrow = grid::arrow(type = "closed", length = grid::unit(6, "mm"), angle = 25), col = "lightgray") +
  geom_path(aes(x, y), col = "lightgray", size = 2) +
  geom_point(aes(x, y, col = col), size = 3) +
  geom_text(aes(ifelse(grepl("dyngen", name), x + 1.5, x), y, label = name, hjust = hjust, vjust = vjust), nudge_x = 0, nudge_y = .1) +
  geom_text(aes(x, y, label = name, colour = col), titdf, size = 6, fontface = "bold") +
  coord_polar() +
  dynplot::theme_graph() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_colour_identity()
ggsave("defence_internal/figures/overview.pdf", g, width = 6, height = 6)


for (i in 2:9) {
  g <- ggplot(df) +
    geom_segment(aes(x = x, xend = x, y = 0, yend = 1.05), titdf, colour = "lightgray", linetype = "dashed") +
    geom_path(aes(x, y, group = gr), arrdf, arrow = grid::arrow(type = "closed", length = grid::unit(6, "mm"), angle = 25), col = "lightgray") +
    geom_path(aes(x, y), col = "lightgray", size = 2) +
    geom_point(aes(x, y, col = col), size = ifelse(grepl(paste0("^§", i), df$name), 10, 3)) +
    coord_polar() +
    dynplot::theme_graph() +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_colour_identity()
  ggsave(paste0("defence_internal/figures/overview_ch", i, ".pdf"), g, width = 2, height = 2)
}

g <- ggplot(df) +
  geom_segment(aes(x = x, xend = x, y = 0, yend = 1.05), titdf, colour = "lightgray", linetype = "dashed") +
  geom_path(aes(x, y, group = gr), arrdf, arrow = grid::arrow(type = "closed", length = grid::unit(6, "mm"), angle = 25), col = "lightgray") +
  geom_path(aes(x, y), col = "lightgray", size = 2) +
  geom_point(aes(x, y, col = col), size = 8) +
  coord_polar() +
  dynplot::theme_graph() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_colour_identity()
ggsave(paste0("defence_internal/figures/overview_all.pdf"), g, width = 2, height = 2)

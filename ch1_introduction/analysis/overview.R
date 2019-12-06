library(tidyverse)
library(scales)

set.seed(1)

df <- tribble(
  ~x, ~y, ~name, ~hjust, ~vjust,
  0, 0, "§2: dyngen", .5, 1,
  .4, 1, "§3: dynbenchmark", 0, .5,
  # .6, 1, "§4: self-assessment\nin trajectory inference", 0, .5,
  .7, 1, "§4: SCORPIUS", 0, .5,
  1, 1, "§5: dyno", 0, .5,
  1.8, 1, "§6: bred", 1, .5,
  2, 1, "§7: incgraph", 1, .5,
  3, 1, "§8: guidelines\nfor benchmarking", .5, 0
)
titdf <- tibble(
  x = 0:2,
  y = c(1.6, 1.8, 1.8),
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
  2, .6, 1,
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
  theme(plot.margin=margin(0,0,0,0)) +
  scale_colour_identity()
ggsave("fig/overview.pdf", g, width = 6, height = 4)



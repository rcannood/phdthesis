library(tidyverse)

extra_evals <-
  tribble(
    ~tool_id, ~num_methods, ~real_datasets, ~simulated_datasets,
    "monocle1", 0, 0, 0,
    "wanderlust", 0, 0, 0,
    "scuba", 0, 0, 0,
    "sincell", 0, 0, 0,
    "nbor", 0, 0, 0,
    "cycler", 0, 0, 0,
    "oscope", 0, 0, 0,
    "waterfall", 0, 0, 0,
    "gpseudotime", 0, 0, 0,
    "embeddr", 0, 0, 0,
    "eclair", 0, 0, 0,
    "dpt", 2, 1, 0,
    "pseudogp", 0, 0, 0,
    "slicer", 0, 0, 0,
    "scell", 0, 0, 0,
    "wishbone", 0, 0, 0,
    "tscan", 5, 3, 0, #sim?
    "scoup", 3, 2, 0, #sim?
    "delorean", 0, 0, 0,
    "raceid_stemid", 0, 0, 0,
    "ouija", 0, 0, 0,
    "mpath", 0, 0, 0,
    "celltree", 3, 1, 0, #sim?
    "wavecrest", 0, 0, 0,
    "stemnet", 0, 0, 0,
    "scimitar", 3, 2, 0, #sim?
    "scorpius", 4, 10, 0,
    "scent", 3, 5, 0,
    "slice", 0, 0, 0,
    "topslam", 0, 0, 0,
    "monocle", 4, 1, 0,
    "gpfates", 0, 0, 0,
    "mfa", 0, 0, 0,
    "tasic", 0, 0, 0,
    "somsc", 0, 0, 0,
    "slingshot", 5, 0, 3,
    "sctda", 4, 0, 2,
    "uncurl", 0, 0, 0,
    "recat", 6, 10, 0,
    "forks", 8, 6, 0,
    "matcher", 1, 1, 1,
    "phenopath", 0, 0, 0,
    "hopland", 7, 6, 5,
    "soptsc", 3, 4, 0,
    "pba", 0, 0, 0,
    "bgp", 3, 0, 1,
    "brgps", 7, 2, 0,
    "wot", 0, 0, 0,
    "treetop", 0, 0, 0,
    "paga", 0, 0, 0,
    "fateid", 0, 0, 0,
    "pcreode", 0, 0, 0,
    "icpsc", 4, 3, 0,
    "grandprix", 2, 1, 0,
    "cshmm", 0, 0, 0,
    "calista", 3, 4, 1,
    "scepath", 5, 4, 0,
    "merlot", 6, 3, 3,
    "gpseudorank", 0, 0, 0,
    "cellrouter", 5, 1, 0,
    "densitypath", 4, 3, 0,
    "topographer", 0, 0, 0,
    "stream", 11, 0, 1,
    "elpigraph", 0, 0, 0,
    "urd", 0, 0, 0,
    "celltrails", 0, 0, 0,
    "ddd", 0, 0, 0,
    "palantir", 0, 0, 0,
    "confess", 0, 0, 0,
    "graphddp", 0, 0, 0,
    "monocle3", 0, 0, 0,
    "psupertime", 4, 5, 0,
    "cyclum", 3, 4, 0,
    "sinova", 0, 0, 0,
    "gpseudoclust", 4, 0, 2,
    "pseudodynamics", 0, 0, 0
  )
tools <-
  read_rds(dynbenchmark::result_file("tools.rds", "03-methods")) %>%
  mutate(
    tool_name = ifelse(tool_id == "monocle", "Monocle 2", ifelse(is.na(tool_name), method_name, tool_name))
  ) %>%
  inner_join(extra_evals, by = c("tool_id")) %>%
  mutate(
    eval_answer = num_methods > 0,
    num_datasets = real_datasets + simulated_datasets,
    peer_reviewed = factor(ifelse(is.na(manuscript_publication_date), "No", "Yes"), levels = c("Yes", "No")),
    manuscript_date2 = ifelse(is.na(manuscript_publication_date), manuscript_preprint_date, manuscript_publication_date) %>% as.Date(origin = "1970-01-01")
  )


pubs <-
  tibble(
    x = as.Date(c("2018-03-05", "2019-04-01"), origin = "1970-01-01"),
    y = 10,
    text = c("bioRxiv", "NBT")
  )
g <- ggplot() +
  geom_point(aes(manuscript_date, num_methods, size = num_datasets, colour = peer_reviewed), tools %>% arrange(desc(num_datasets)), shape = 21) +
  geom_vline(aes(xintercept = x), pubs, linetype = "dashed", colour = "darkgray") +
  geom_text(aes(x, 10, label = text), pubs, colour = "darkgray", hjust = 1, nudge_x = -15) +
  theme_classic() +
  labs(x = "Date", y = "Number of methods", size = "Number of datasets", colour = "Peer reviewed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  scale_colour_brewer(palette = "Set1") +
  scale_size_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 12))

g <- ggplot() +
  geom_point(aes(manuscript_date, num_datasets, size = num_methods, colour = peer_reviewed), tools %>% arrange(desc(num_methods)), shape = 21) +
  geom_vline(aes(xintercept = x), pubs, linetype = "dashed", colour = "darkgray") +
  geom_text(aes(x, 8, label = text), pubs, colour = "darkgray", hjust = 1, nudge_x = -15) +
  theme_classic() +
  labs(x = "Date", y = "Number of datasets", size = "Number of methods", colour = "Peer reviewed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  scale_colour_brewer(palette = "Set1") +
  scale_size_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 12))

ggsave("fig/benchmarks.pdf", g, width = 6, height = 3)
# ggplot(tools, aes(manuscript_date, num_methods)) +
#   geom_point() +
#   geom_vline(aes(xintercept = x), tibble(x = as.Date(c("2019-04-01", "2018-03-05"))), linetype = "dashed", colour = "darkgray") +
#   theme_classic() +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y")
# ggplot(tools, aes(manuscript_date, num_datasets)) +
#   geom_point() +
#   geom_vline(aes(xintercept = x), tibble(x = as.Date(c("2019-04-01", "2018-03-05"))), linetype = "dashed", colour = "darkgray") +
#   theme_classic() +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y")
# ggplot(tools, aes(manuscript_date, num_datasets)) +
#   geom_text(aes(label = num_methods)) +
#   geom_vline(aes(xintercept = x), tibble(x = as.Date(c("2019-04-01", "2018-03-05"))), linetype = "dashed", colour = "darkgray") +
#   theme_classic() +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y")
#
# ggplot(tools, aes(num_datasets, num_methods)) +
#   geom_point() +
#   ggrepel::geom_text_repel(aes(label = tool_name), tools %>% filter(num_methods > 0)) +
#   theme_classic()
#
# tools_dat <-
#   tools %>%
#   select(tool_id, tool_name, manuscript_date, num_methods, num_datasets) %>%
#   gather(variable, value, num_methods, num_datasets) %>%
#   mutate(
#     date_quart = lubridate::floor_date(manuscript_date, "quarter") + 30,
#     group = case_when(
#       value == 0 ~ "0",
#       value == 1 ~ "1",
#       value <= 3 ~ "2-3",
#       value <= 5 ~ "4-5",
#       value <= 8 ~ "6-8",
#       TRUE ~ "9-11"
#     ) %>% factor(levels = c("0", "1", "2-3", "4-5", "6-8", "9-11"))
#   )
# ggplot(tools_dat) +
#   geom_histogram(aes(date_quart, fill = forcats::fct_rev(group)), binwidth = 365) +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
#   scale_fill_brewer(palette = "Blues", direction = -1) +
#   theme_classic() +
#   facet_wrap(~variable) +
#   guides(fill = guide_legend(reverse = TRUE))
#
# summ_tools <- tools %>%
#   select(num_methods, num_datasets) %>%
#   gather(variable, value) %>%
#   mutate(value = case_when(value == 0 ~ "0", value == 1 ~ "1", value <= 3 ~ "2-3", value <= 5 ~ "4-5", value <= 8 ~ "6-8", TRUE ~ "9-11") %>% factor(levels = c("0", "1", "2-3", "4-5", "6-8", "9-11"))) %>%
#   group_by(variable, value) %>%
#   summarise(pct = n()) %>%
#   mutate(pct = pct / sum(pct)) %>%
#   ungroup()
# ggplot(summ_tools) +
#   geom_bar(aes(variable, pct, fill = value, group = forcats::fct_rev(value)), stat = "identity") +
#   scale_y_continuous(breaks = c(0, .25, .5, .75, 1), labels = rev(c("100%", "75%", "50%", "25%", "0%"))) +
#   # scale_fill_brewer(palette = "Paired") +
#   scale_fill_brewer(palette = "Blues") +
#   # scale_fill_gradientn(colours = RColorBrewer::brewer.pal(8, "Blues")[-1:-2]) +
#   coord_flip() +
#   theme_classic()
#
#
# ggplot(tools, aes(manuscript_date, log10(manuscript_citations + 1))) +
#   geom_point(aes(colour = num_methods, size = real_datasets + simulated_datasets)) +
#   ggrepel::geom_text_repel(aes(label = tool_name), colour = "darkgray") +
#   theme_bw() +
#   viridis::scale_color_viridis(direction = -1) +
#   scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
#   labs(size = "# Datasets", colour = "# Methods", shape = "Benchmarked", x = "Time", y = "Log Citations + 1")
#
#
#
#
# ggplot(tools, aes(manuscript_date, log10(manuscript_citations + 1))) +
#   geom_point(aes(size = num_methods, colour = real_datasets + simulated_datasets)) +
#   geom_text(aes(label = tool_name), colour = "darkgray", nudge_y = .1) +
#   theme_bw() +
#   viridis::scale_color_viridis(direction = -1) +
#   scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
#   labs(colour = "# Datasets", size = "# Methods", shape = "Benchmarked", x = "Time", y = "Log Citations + 1")
#
# ggplot(tools, aes(manuscript_date, num_methods, colour = factor(eval_answer))) +
#   geom_point(aes()) +
#   geom_text(aes(label = tool_id), nudge_y = .1) +
#   theme_bw()
# ggplot(tools, aes(manuscript_date, real_datasets, colour = factor(eval_answer))) +
#   geom_point(aes()) +
#   geom_text(aes(label = tool_id), nudge_y = .1) +
#   theme_bw()
# ggplot(tools, aes(manuscript_date, simulated_datasets, colour = factor(eval_answer))) +
#   geom_point(aes()) +
#   geom_text(aes(label = tool_id), nudge_y = .1) +
#   theme_bw()

dates <- seq(as.Date("2014-01-01"), Sys.Date(), by = 1)
tools_cum <-
  tools %>%
  select(tool_id, tool_name, manuscript_date, num_methods, num_datasets) %>%
  gather(variable, value, num_methods, num_datasets) %>%
  mutate(
    group = case_when(
      value == 0 ~ "0",
      value == 1 ~ "1",
      value <= 3 ~ "2-3",
      value <= 5 ~ "4-5",
      value <= 8 ~ "6-8",
      TRUE ~ "9-11"
    ) %>% factor(levels = c("0", "1", "2-3", "4-5", "6-8", "9-11"))
  ) %>%
  arrange(manuscript_date) %>%
  crossing(groupings = c("0", "1", "2-3", "4-5", "6-8", "9-11")) %>%
  group_by(variable, groupings) %>%
  do({
    df <- .
    df <- df %>% filter(group == groupings)
    if (nrow(df) > 0) {
      tib <- tibble(
        variable = df$variable[[1]],
        groupings = df$groupings[[1]],
        manuscript_date = dates,
        count = map_int(manuscript_date, function(da) sum(df$manuscript_date == da)),
        cum = cumsum(count)
      )
    } else {
      NULL
    }
  }) %>%
  ungroup() %>%
  mutate(facet = c(num_datasets = "Number of datasets", num_methods = "Number of methods")[variable])

g1 <- ggplot(tools_cum %>% filter(variable == "num_datasets")) +
  geom_area(aes(x = manuscript_date, y = cum, fill = forcats::fct_rev(groupings))) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "# Articles", fill = "# Datasets used\nin self-benchmark", tag = "A") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Blues", direction = -1)
g2 <- ggplot(tools_cum %>% filter(variable == "num_methods")) +
  geom_area(aes(x = manuscript_date, y = cum, fill = forcats::fct_rev(groupings))) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "# Articles", fill = "# Methods compared\nin self-benchmark", tag = "B") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Blues", direction = -1)
ggsave("fig/self_assessment.pdf", patchwork::wrap_plots(g1, g2, ncol = 1), width = 6, height = 6)

^tools_pie <-
  tools %>%
  select(manuscript_citations, eval_answer) %>%
  arrange(eval_answer, manuscript_citations) %>%
  mutate(name = forcats::fct_inorder(paste0("row_", row_number()))) %>%
  group_by(eval_answer) %>%
  mutate(
    colour = colorRampPalette(RColorBrewer::brewer.pal(8, ifelse(eval_answer[[1]], "Blues", "Reds")))(1000)[round(cumsum(manuscript_citations) / sum(manuscript_citations) * 999) + 1]
  ) %>%
  ungroup()
ggplot(tools_pie) +
  geom_bar(aes(1, manuscript_citations, fill = colour, group = name), stat = "identity", colour = "black") +
  coord_polar("y") +
  scale_fill_identity()
g <- ggplot(tools_pie) +
  geom_rect(aes(xmin = 0, xmax = 0, ymin = 0, ymax = 0, fill = group), tibble(group = c("Self-assessment", "No self-assessment"))) +
  scale_fill_manual(values = c("Self-assessment" = "#084594", "No self-assessment" = "#99000D")) +
  labs(fill = "Group", x = NULL, y = "Number of citations") +
  ggnewscale::new_scale_fill() +
  geom_bar(aes(0, manuscript_citations, fill = colour, group = name), stat = "identity") +
  coord_polar("y") +
  scale_fill_identity() +
  theme_classic() +
  scale_x_continuous(breaks = NULL)#+
  # scale_y_continuous(breaks = 1:4 / 4 * 5500)
g
ggsave("fig/citations.pdf", g, width = 5, height = 3)

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

# cumulative plot
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
  labs(x = NULL, y = "# Articles", fill = "# Datasets used\nin self-assessment", tag = "A") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Blues", direction = -1)
g2 <- ggplot(tools_cum %>% filter(variable == "num_methods")) +
  geom_area(aes(x = manuscript_date, y = cum, fill = forcats::fct_rev(groupings))) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "# Articles", fill = "# Methods compared\nin self-assessment", tag = "B") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Blues", direction = -1)
ggsave("fig/self_assessment.pdf", patchwork::wrap_plots(g1, g2, ncol = 1), width = 6, height = 6)

# citation pie
tools_pie <-
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
  scale_x_continuous(breaks = NULL)
g
ggsave("fig/citations.pdf", g, width = 5, height = 3)


datasets <- dynbenchmark::load_datasets()
dynbenchmark::list_datasets("real/gold/germlin")
real_datasets <- read_rds(dynbenchmark::result_file("metadata.rds", "01-datasets/01-real"))

dates <- seq(as.Date("2014-01-01"), as.Date("2018-12-31"), by = 1)
trajtyps <- dynwrap::trajectory_types %>%
  filter(!id %in% c("convergence", "acyclic_graph")) %>%
  mutate(
    label = dynbenchmark::label_long(id)
  )
datasets_cum <-
  real_datasets %>%
  filter(trajectory_type %in% trajtyps$id) %>%
  mutate(trajectory_type = factor(trajectory_type, levels = trajtyps$id)) %>%
  group_by(trajectory_type) %>%
  do({
    df <- .
    if (nrow(df) > 0) {
      tibble(
        trajectory_type = df$trajectory_type[[1]],
        date = dates,
        count = map_int(date, function(da) sum(df$date == da)),
        cum = cumsum(count)
      )
    } else {
      NULL
    }
  }) %>%
  ungroup() %>%
  mutate(
    trajectory_type_label = factor(dynbenchmark::label_long(trajectory_type), levels = trajtyps$label)
  )
g <- ggplot(datasets_cum) +
  geom_area(aes(x = date, y = cum, fill = forcats::fct_rev(trajectory_type_label))) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "# Datasets", fill = "Trajectory type") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2019-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 120), expand = c(0, 0)) +
  scale_fill_manual(values = trajtyps %>% select(label, colour) %>% deframe())
ggsave("fig/datasets.pdf", g, width = 6, height = 2.5)


bench_out <- read_rds(dynbenchmark::result_file("benchmark_results_normalised.rds", "06-benchmark"))
#boda <- bench_out$data_aggregations %>%
#  filter(dataset_trajectory_type == "overall", dataset_source != "mean") %>%
boda <- bench_out$data %>%
  filter(method_id %in% dynbenchmark::load_methods()$method_id) %>%
  # filter(dataset_source %in% c("real/gold", "real/silver", "synthetic/dyntoy", "synthetic/dyngen")) %>%
  filter(dataset_source %in% c("real/gold", "real/silver", "synthetic/dyngen")) %>%
  mutate(dataset_source2 = gsub("/.*", "", dataset_source)) %>%
  group_by(dataset_source2, method_id, method_name) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup()
ggplot(boda %>% select(dataset_source2, overall, method_id, method_name) %>% spread(dataset_source2, overall)) +
  geom_point(aes(real, synthetic)) +
  geom_abline(intercept = 0, slope = 1) +
  expand_limits(x = c(0, 1), y = c(0, 1))

boda <- bench_out$data %>%
  filter(
    method_id %in% dynbenchmark::load_methods()$method_id
  ) %>%
  mutate(dataset_source2 = gsub("/.*", "", dataset_source)) %>%
  group_by(dataset_source2, dataset_trajectory_type, method_id, method_name) %>%
  summarise(overall = mean(overall), n = n()) %>%
  ungroup() %>%
  filter(n > 5, !dataset_trajectory_type %in% c("acyclic_graph", "convergence")) %>%
  select(-n) %>%
  spread(dataset_source2, overall) %>%
  na.omit()
ggplot(boda) +
  geom_point(aes(real, synthetic, colour = dataset_trajectory_type)) +
  geom_abline(intercept = 0, slope = 1) +
  expand_limits(x = c(0, 1), y = c(0, 1)) +
  coord_equal()


ggplot(boda) +
  geom_density_2d(aes(real, synthetic)) +
  geom_abline(intercept = 0, slope = 1) +
  expand_limits(x = c(0, 1), y = c(0, 1))



boda <- bench_out$data %>%
  filter(
    method_id %in% dynbenchmark::load_methods()$method_id
  ) %>%
  mutate(dataset_source2 = ifelse(grepl("real", dataset_source), "real", dataset_source)) %>%
  group_by(dataset_source2, dataset_trajectory_type, method_id, method_name) %>%
  summarise(overall = mean(overall), n = n()) %>%
  ungroup() %>%
  filter(n > 5, !dataset_trajectory_type %in% c("acyclic_graph", "convergence")) %>%
  select(-n) %>%
  spread(dataset_source2, overall) %>%
  gather(synthetic_source, synthetic_score, -dataset_trajectory_type:-real)
ggplot(boda) +
  geom_point(aes(real, synthetic_score, colour = dataset_trajectory_type)) +
  facet_wrap(~synthetic_source) +
  geom_abline(intercept = 0, slope = 1) +
  expand_limits(x = c(0, 1), y = c(0, 1))

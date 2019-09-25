library(tidyverse)


none <- character(0)


extra_evals <-
  tribble(
    ~tool_id, ~methods, ~real_datasets, ~simulated_datasets, ~metrics,
    "monocle1", none, none, none, none,
    "wanderlust", none, none, none, none,
    "scuba", none, none, none, none,
    "sincell", none, none, none, none,
    "nbor", none, none, none, none,
    "cycler", none, none, none, none,
    "oscope", none, none, none, none,
    "waterfall", none, none, none, none,
    "gpseudotime", none, none, none, none,
    "embeddr", none, none, none, none,
    "eclair", none, none, none, none,
    "dpt", c("dpt", "wishbone", "monocle1"), c("moignard", "klein", "paul"), c("synthetic_grn"), c("pseudotime_correlation", "robustness"),
    "pseudogp", none, none, none, none,
    "slicer", none, none, none, none,
    "scell", none, none, none, none,
    "wishbone", none, none, none, none,
    "tscan", c("monocle1", "tscan", "waterfall", "scuba", "wanderlust"), c("trapnell", "amit", "shin"), none, "pseudotime_pos",
    "scoup", c("scoup", "monocle1", "tscan"), c("kouno", "moignard", "shalek"), none, "pseudotime_pis",
    "delorean", none, none, none, none,
    "raceid_stemid", none, none, none, none,
    "ouija", none, none, none, none,
    "mpath", none, none, none, none,
    "celltree", c("monocle1", "tscan", "celltree"), "trapnell", none, "pseudotime_unknown",
    "wavecrest", none, none, none, none,
    "stemnet", none, none, none, none,
    "scimitar", c("scimitar", "monocle1", "wanderlust"), none, c("synthetic_custom", "synthetic_custom"), "pseudotime_correlation",
    "scorpius", c("scorpius", "wanderlust", "monocle1", "waterfall"), c("schlitzer", "buettner", "shalek", "trapnell", "kowalczyk"), none, c("pseudotime_cos", "robustness_cva"),
    "scent", c("scent", "slice", "stemid"), c("chu", "trapnell", "treutlein"), none, c("pseudotime_wilcox", "pseudotime_auc"),
    "slice", none, none, none, none,
    "topslam", c("monocle1", "wishbone", "topslam"), none, "synthetic_zwiessele", "pseudotime_correlation",
    "monocle", c("monocle1", "monocle", "dpt", "wishbone"), "paul", none, c("pseudotime_correlation", "branch_ari"),
    "gpfates", none, none, none, none,
    "mfa", none, none, none, none,
    "tasic", none, none, none, none,
    "somsc", none, none, none, none,
    "slingshot", c("slingshot", "monocle1", "monocle", "dpt", "tscan"), none, paste0("splatter", 1:5), c("treepseudotime_correlation"),
    "sctda", c("sctda", "wishbone", "slicer", "dpt"), none, "synthetic_rizvi", "pseudotime_correlation",
    "uncurl", none, none, none, none,
    "recat", c("recat", "scuba", "monocle1", "tscan", "wishbone", "dpt"), "buettner", none, c("pseudotime_correlation", "pseudotime_custom"),
    "forks", c("forks", "monocle", "scuba", "tscan", "waterfall", "dpt", "gpfates", "slicer"), c("windram", "deng", "guo", "klein", "amit", "petropoulos"), none, c("pseudotime_correlation", "robustness_stdev"),
    "matcher", "matcher", "angelmueller", "synthetic_custom", "pseudotime_correlation",
    "phenopath", none, none, none, none,
    "hopland", c("hopland", "wanderlust", "monocle1", "topslam", "scuba", "wishbone", "dpt"), c("guo", "deng", "yan", "amit", "islam"), "synthetic_zwiessele", "pseudotime_correlation",
    "soptsc", c("soptsc", "monocle", "dpt"), c("guo", "klein", "shalek"), none, "pseudotime_correlation",
    "pba", none, none, none, none,
    "brgps", c("brgps", "grandprix", "monocle", "scuba", "slicer", "tscan", "wishbone"), c("guo", "guo"), none, "pseudotime_correlation",
    "wot", none, none, none, none,
    "treetop", none, none, none, none,
    "paga", none, none, none, none,
    "fateid", none, none, none, none,
    "pcreode", none, none, none, none,
    "icpsc", c("icpsc", "wishbone", "monocle", "dpt"), c("sun", "trapnell", "yao"), none, "pseudotime_correlation",
    "grandprix", c("delorean", "grandprix"), c("windram"), none, "pseudotime_correlation",
    "cshmm", none, none, none, none,
    "calista", c("monocle", "calista", "dpt"), c("moignard", "bargaje", "treutlein", "chu"), c("synthetic_custom"), "pseudotime_correlation",
    "scepath", c("scepath", "monocle1", "monocle", "tscan", "dpt"), c("yan", "treutlein", "trapnell"), none, c("pseudotime_correlation", "robustness_correlation"),
    "merlot", c("merlot", "dpt", "slicer", "monocle", "slingshot", "tscan"), c("paul", "guo", "velten"), c("prosstt1", "prosstt2", "splatter"), c("branch_mi", "pseudotime_correlation"),
    "gpseudorank", none, none, none, none,
    "cellrouter", c("monocle", "dpt", "wishbone", "waterfall"), c("paul", "olsson"), none, "internal_autocorrelation",
    "densitypath", c("monocle", "wishbone", "dpt", "densitypath"), c("petropoulos"), c("synthetic_moon", "synthetic_zwiessele"), c("branch_ari", "pseudotime_correlation"),
    "topographer", none, none, none, none,
    "stream", c("stream", "sctda", "wishbone", "slicer", "monocle", "dpt", "tscan", "scuba", "mpath", "gpfates"), none, "synthetic_rizvi", "pseudotime_correlation",
    "elpigraph", none, none, none, none,
    "urd", none, none, none, none,
    "celltrails", none, none, none, none,
    "ddd", none, none, none, none,
    "palantir", none, none, none, none,
    "confess", none, none, none, none,
    "graphddp", none, none, none, none,
    "monocle3", none, none, none, none,
    "psupertime", c("monocle", "slingshot", "psupertime"), c("enge", "qiu", "petropoulos", "li", "treutlein"), none, "pseudotime_correlation",
    "cyclum", c("cyclum", "recat"), c("buettner", "mcdavid", "mcdavid", "mcdavid"), none, c("cluster_accuracy"),
    "sinova", none, none, none, none,
    "gpseudoclust", c("gpseudoclust", "monocle", "delorian", "slicer"), c("sasagawa", "shalek"), c("synthetic_strauss", "synthetic_strauss"), c("cluster_ari", "cluster_fmi", "cluster_nmi"),
    "pseudodynamics", none, none, none, none
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

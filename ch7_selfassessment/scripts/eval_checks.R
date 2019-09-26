library(tidyverse)
library(googlesheets)
library(rlang)


# COLLECT DATA ------------------------------------------------------------

# If it's your first time running this script, run this:
# gs_auth()
sheet <- gs_key("1Mug0yz8BebzWt8cmEW306ie645SBh_tDHwjVw4OFhlE")

case_study_data <- gs_read(sheet, ws = "case_study") %>% select(tool_id = id, methods, datasets, metrics)

tools <-
  read_rds(dynbenchmark::result_file("tools.rds", "03-methods")) %>%
  mutate(
    tool_id = ifelse(tool_id == "monocle", "monocle2", tool_id),
    tool_name = ifelse(tool_id == "monocle2", "Monocle 2", ifelse(is.na(tool_name), method_name, tool_name)),
    manuscript_doi = ifelse(tool_id == "elpigraph", "arXiv:1804.0758", manuscript_doi)
  ) %>%
  inner_join(case_study_data, by = c("tool_id")) %>%
  mutate_at(c("methods", "datasets", "metrics"), ~ strsplit(. %|% "", ",")) %>%
  mutate(
    synthetic_datasets = map(datasets, keep, ~ grepl("synthetic_", .)),
    real_datasets = map(datasets, discard, ~ grepl("synthetic_", .)),
    num_methods = map_int(methods, length),
    num_datasets = map_int(datasets, length),
    num_metrics = map_int(metrics, length),
    self_assessment = num_methods > 0,
    peer_reviewed = !is.na(manuscript_publication_date),
    has_preprint = !is.na(manuscript_preprint_date),
    manuscript_date = ifelse(has_preprint, manuscript_preprint_date, manuscript_publication_date) %>% as.Date(origin = "1970-01-01"),
    manuscript_date2 = ifelse(peer_reviewed, manuscript_publication_date, manuscript_preprint_date) %>% as.Date(origin = "1970-01-01"),
    metric_group = map_chr(metrics, function(w) {
      x <- gsub("_.*", "", w)
      case_when(
        length(x) == 0 ~ "None",
        "branch" %in% x && "pseudotime" %in% x ~ "Cluster & pseudotime",
        "cluster" %in% x ~ "Cluster",
        any(c("pseudotime", "internal", "pseudotime") %in% x) ~ "Pseudotime",
        "cyclepseudotime" %in% x ~ "Pseudotime for cycles",
        "treepseudotime" %in% x ~ "Pseudotime for trees"
      )
    }) %>% factor(levels = c("None", "Pseudotime", "Cluster", "Cluster & pseudotime", "Pseudotime for cycles", "Pseudotime for trees"))
  )

if (save_tsv) {
  # # Only run this if you are prepared to clean up the mess if something goes wrong ;)
  # overwrite <-
  #   tools %>%
  #   arrange(manuscript_date) %>%
  #   select(
  #     id = tool_id,
  #     name = tool_name,
  #     preprint_date = manuscript_preprint_date,
  #     pub_date = manuscript_publication_date,
  #     doi = manuscript_doi,
  #     google_scholar_id = manuscript_google_scholar_cluster_id,
  #     methods,
  #     metrics,
  #     datasets
  #   ) %>%
  #   mutate_at(c("methods", "metrics", "datasets"), ~ map_chr(., paste, collapse = ",")) %>%
  #   mutate_at(c("preprint_date", "pub_date"), function(x) ifelse(is.na(x), "", as.character(x)))
  # gs_edit_cells(sheet, ws = "case_study", input = overwrite, anchor = "A1", col_names = TRUE, trim = TRUE)
}



# INTRODUCTION ------------------------------------------------------------

tools
tools$self_assessment %>% mean
tools %>% filter(!is.na(manuscript_publication_date)) %>% pull(self_assessment) %>% mean
tools %>% filter(!is.na(manuscript_preprint_date)) %>% pull(self_assessment) %>% mean
tools %>% filter(num_methods > 5, num_datasets > 5)

tools %>% filter(num_datasets > 4)

tools %>% filter(map_lgl(synthetic_datasets, ~length(.) > 0)) %>% transmute(real = map_lgl(real_datasets, ~length(.) > 0)) %>% table
tools %>% filter(map_lgl(synthetic_datasets, ~length(.) > 0)) %>% select(tool_id, synthetic_datasets) %>% unnest(syn = synthetic_datasets) %>% as.data.frame

tools %>% select(metrics, tool_id) %>% unnest(metric = metrics) %>% mutate(metric = gsub("_.*", "", metric)) %>% unique() %>% group_by(metric) %>% summarise(n = n(), m = paste(tool_id, collapse = ", "))

# cumulative plot
dates <- seq(as.Date("2014-01-01"), Sys.Date(), by = 1)


# pct cumulatives
pct_cum <- bind_rows(
  tools %>% mutate(group = "All", date = manuscript_date),
  tools %>% filter(has_preprint) %>% mutate(group = "Pre-print", date = manuscript_preprint_date),
  tools %>% filter(peer_reviewed) %>% mutate(group = "Peer-reviewed", date = manuscript_publication_date)
) %>%
  group_by(group) %>%
  do({
    df <- .
    df2 <- df %>% filter(self_assessment)
    tibble(
      group = df$group[[1]],
      date = dates,
      num_methods = cumsum(map_int(date, function(da) sum(df$date == da))),
      num_evals = cumsum(map_int(date, function(da) sum(df2$date == da))),
      pct_eval = num_evals / num_methods
    )
  }) %>%
  ungroup()

g1 <- ggplot(pct_cum) +
  geom_step(aes(x = date, y = num_methods, colour = group)) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "# Articles", colour = "Group", tag = "A") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_colour_brewer(palette = "Set1")
g2 <- ggplot(pct_cum) +
  geom_step(aes(x = date, y = pct_eval, colour = group)) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "% self-assessment", colour = "Group", tag = "B") +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0)) +
  scale_colour_brewer(palette = "Set1")

g3 <- ggplot(tools) +
  geom_point(aes(manuscript_date, num_datasets)) +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_bw() +
  labs(x = NULL, y = "# Datasets", colour = "Metric group", tag = "C")

g4 <- ggplot(tools) +
  geom_point(aes(manuscript_date, num_methods)) +
  scale_x_date(limits = as.Date(c("2014-01-01", "2020-01-01")), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_bw() +
  labs(x = NULL, y = "# Methods", colour = "Metric group", tag = "D")

g <- patchwork::wrap_plots(g1, g2 + theme(legend.position = "none"), g3, g4 + theme(legend.position = "none"), ncol = 1)
ggsave("fig/self_assessment.pdf", g, width = 6, height = 6)


# DATASETS ----------------------------------------------------------------
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


# METRICS -----------------------------------------------------------------
metric_summ <-
  tools %>%
  mutate(metric_group = forcats::fct_infreq(metric_group)) %>%
  group_by(metric_group) %>%
  summarise(n = n()) %>%
  mutate(
    xmax = cumsum(n),
    xmin = xmax - n,
    xmid = (xmin + xmax) / 2
  )
g <- ggplot(metric_summ) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = metric_group)) +
  geom_text(aes(xmid, 1.1, label = n)) +
  coord_polar() +
  theme_classic() +
  dynplot::theme_graph() +
  scale_fill_manual(values = c("lightgray", RColorBrewer::brewer.pal(5, "Set3"))) +
  labs(fill = "Metric group")
g
ggsave("fig/metrics.pdf", g, width = 5, height = 3)

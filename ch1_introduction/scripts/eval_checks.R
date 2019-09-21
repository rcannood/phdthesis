library(tidyverse)

# qc_checks <- read_rds(dynbenchmark::result_file("qc_checks.rds", "03-methods")) %>%
#   filter(aspect_id == "evaluation", check_id == 48)
evals <- read_rds(dynbenchmark::result_file("tool_qc.rds", "03-methods")) %>%
  filter(check_id == 48) %>%
  select(tool_id, eval_answer = answer)
tools %>% arrange(manuscript_date) %>% pull(tool_id)
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
    "sinova", 0, 0, 0
  )
tools <-
  read_rds(dynbenchmark::result_file("tools.rds", "03-methods")) %>%
  mutate(
    tool_name = ifelse(tool_id == "monocle", "Monocle 2", ifelse(is.na(tool_name), method_name, tool_name))
  ) %>%
  inner_join(extra_evals, by = c("tool_id")) %>%
  mutate(
    eval_answer = num_methods > 0,
    num_datasets = real_datasets + simulated_datasets
  )

ggplot(tools, aes(num_datasets, num_methods)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = tool_name), tools %>% filter(num_methods > 0)) +
  theme_classic()

summ_tools <- tools %>%
  select(num_methods, num_datasets) %>%
  gather(variable, value) %>%
  mutate(value = case_when(value == 0 ~ "0", value == 1 ~ "1", value <= 3 ~ "2-3", value <= 5 ~ "4-5", value <= 8 ~ "6-8", TRUE ~ "9-11") %>% factor(levels = c("0", "1", "2-3", "4-5", "6-8", "9-11"))) %>%
  group_by(variable, value) %>%
  summarise(pct = n()) %>%
  mutate(pct = pct / sum(pct)) %>%
  ungroup()
ggplot(summ_tools) +
  geom_bar(aes(variable, pct, fill = value, group = forcats::fct_rev(value)), stat = "identity") +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), labels = rev(c("100%", "75%", "50%", "25%", "0%"))) +
  # scale_fill_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Blues") +
  # scale_fill_gradientn(colours = RColorBrewer::brewer.pal(8, "Blues")[-1:-2]) +
  coord_flip() +
  theme_classic()


ggplot(tools, aes(manuscript_date, log10(manuscript_citations + 1))) +
  geom_point(aes(colour = num_methods, size = real_datasets + simulated_datasets)) +
  ggrepel::geom_text_repel(aes(label = tool_name), colour = "darkgray") +
  theme_bw() +
  viridis::scale_color_viridis(direction = -1) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
  labs(size = "# Datasets", colour = "# Methods", shape = "Benchmarked", x = "Time", y = "Log Citations + 1")




ggplot(tools, aes(manuscript_date, log10(manuscript_citations + 1))) +
  geom_point(aes(size = num_methods, colour = real_datasets + simulated_datasets)) +
  geom_text(aes(label = tool_name), colour = "darkgray", nudge_y = .1) +
  theme_bw() +
  viridis::scale_color_viridis(direction = -1) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
  labs(colour = "# Datasets", size = "# Methods", shape = "Benchmarked", x = "Time", y = "Log Citations + 1")

ggplot(tools, aes(manuscript_date, num_methods, colour = factor(eval_answer))) +
  geom_point(aes()) +
  geom_text(aes(label = tool_id), nudge_y = .1) +
  theme_bw()
ggplot(tools, aes(manuscript_date, real_datasets, colour = factor(eval_answer))) +
  geom_point(aes()) +
  geom_text(aes(label = tool_id), nudge_y = .1) +
  theme_bw()
ggplot(tools, aes(manuscript_date, simulated_datasets, colour = factor(eval_answer))) +
  geom_point(aes()) +
  geom_text(aes(label = tool_id), nudge_y = .1) +
  theme_bw()

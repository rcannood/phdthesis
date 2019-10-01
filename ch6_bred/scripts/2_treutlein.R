library(tidyverse)
library(dyno)
options(dynwrap_backend = 'r_wrapper')

# download and preprocess dataset
file <- "derived_files/GSE67310.txt.gz"
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67310&format=file&file=GSE67310%5FiN%5Fdata%5Flog2FPKM%5Fannotated%2Etxt%2Egz"
if (!file.exists(file)) {
  download.file(url, file)
}
data <- read_tsv(file, col_types = cols(
  cell_name = "c", assignment = "c", experiment = "c", time_point = "c", .default = "d"
))

cell_info <-
  data %>%
  select(cell_name:time_point) %>%
  rename(cell_id = cell_name)
grouping <-
  cell_info %>%
  select(cell_id, assignment) %>%
  deframe()
expression <-
  data %>%
  select(-cell_name:-time_point) %>%
  as.matrix() %>%
  magrittr::set_rownames(data$cell_name)

alias2eg <- as.list(org.Mm.eg.db::org.Mm.egALIAS2EG) %>% unlist()
feature_info <-
  tibble(
    feature_id = colnames(expression),
    entrez = ifelse(feature_id %in% names(alias2eg), alias2eg[feature_id], NA_character_)
  )

counts <- floor(2^expression-1)
dataset <-
  wrap_expression(
    id = "treutlein_2016",
    counts = counts,
    expression = expression,
    cell_info = cell_info,
    feature_info = feature_info
  ) %>%
  add_grouping(
    grouping = grouping
  )

# run slingshot and postprocess trajectory
model <- infer_trajectory(dataset, ti_slingshot())
plot_dimred(model, label_milestones = TRUE, grouping = dataset$grouping)
model2 <- model %>%
  simplify_trajectory() %>%
  label_milestones(c("3" = "MEF", "2" = "Induced", "5" = "Myocyte", "4" = "Neuron")) %>%
  add_root(root_milestone_id = "3")
model2$dimred_segment_points <-
  model2$dimred_segment_progressions <-
  NULL

# calculate feature importance scores
fimp <- dynfeature::calculate_overall_feature_importance(model2, dataset)

# create visualisations
g <- plot_dimred(model2, label_milestones = TRUE, grouping = dataset$grouping)
ggsave("derived_files/plot_dimred.pdf", g, width = 6, height = 4)

g <- plot_heatmap(model2, dataset, fimp$feature_id[1:150])
ggsave("derived_files/plot_heatmap.pdf", g, width = 12, height = 12)

# select targets, regulators and samples
targets <- fimp %>% filter(importance > .001) %>% pull(feature_id)

genesets <- read_rds("derived_files/data_genesets.rds")
reg_entrezs <- genesets %>% filter(grepl("transcription factor", description)) %>% pull(entrezs) %>% unlist() %>% unique()
regulators <- dataset$feature_info %>% filter(entrez %in% reg_entrezs, feature_id %in% targets) %>% pull(feature_id)

samples <- rownames(expression)

# write data to file
write_rds(lst(dataset, model = model2, fimp, targets, regulators, samples), "derived_files/data.rds", compress = "gz")

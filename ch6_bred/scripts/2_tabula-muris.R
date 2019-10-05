library(tidyverse)

dest_dir <- "derived_files/tabula_muris/"
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

data_rnaseq <- paste0(dest_dir, "data_rnaseq/")
data_rds <- paste0(dest_dir, "data_rds/")

# zip_file <- paste0(dest_dir, "data_rnaseq.zip")
# download.file("https://ndownloader.figshare.com/articles/5968960/versions/3", zip_file)
# unzip(zip_file, exdir = data_rnaseq)
# unzip(paste0(data_rnaseq, "droplet.zip"), exdir = data_rnaseq)

# cell_info <- read_csv(
#   paste0(data_rnaseq, "annotations_droplet.csv"),
#   col_types = cols(
#     .default = "c",
#     tissue_tSNE_1 = "d",
#     tissue_tSNE_2 = "d",
#     subsetA = "l",
#     subsetB = "l",
#     subsetC = "l",
#     subsetD = "l"
#   )
# )
#
# write_rds(cell_info, paste0(data_rds, "cell_info.rds"), compress = "gz")
#
# droplet <- paste0(data_rnaseq, "/droplet")
# samples <- list.files(droplet, full.names = TRUE)
# counts <- pbapply::pblapply(samples, function(folder_10x) {
#   channel <- basename(folder_10x) %>% gsub(".*-", "", .)
#
#   counts <- Seurat::Read10X(folder_10x) %>% Matrix::t()
#   rownames(counts) <- paste0(channel, "_", rownames(counts))
#
#   counts <- counts[rownames(counts) %in% cell_info$cell, ]
#
#   counts
# }) %>% do.call(rbind, .)
# counts <- counts[cell_info$cell, ] # reorder cells
#
# write_rds(counts, paste0(data_rds, "counts.rds"), compress = "gz")
#
# tpm <- counts / Matrix::rowSums(counts) * 1000000
# # write_rds(tpm, paste0(data_rds, "tpm.rds"), compress = "gz")
# write_rds(tpm, paste0(data_rds, "tpm.rds"))
#
# expression <- tpm
# expression@x <- log2(expression@x + 1)
# # write_rds(expression, paste0(data_rds, "expression.rds"), compress = "gz")
# write_rds(expression, paste0(data_rds, "expression.rds"))


# FILTER DATA -------------------------------------------------------------
data_rds <- "derived_files/tabula_muris/data_rds/"

if (!file.exists(paste0(data_rds, "data.rds"))) {
  counts <- read_rds(paste0(data_rds, "counts.rds"))
  # expression <- read_rds(paste0(data_rds, "expression.rds"))
  cell_info <- read_rds(paste0(data_rds, "cell_info.rds"))

  seu <-
    Seurat::CreateSeuratObject(
      Matrix::t(counts),
      min.cells = 10,
      min.features = 500,
      meta.data = cell_info %>% column_to_rownames("cell")
    ) %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()
  expression <- Matrix::t(seu@assays$RNA@data)
  feature_info <- seu@assays$RNA@meta.features %>% rownames_to_column("feature_id") %>% as_tibble()
  sample_info <- seu@meta.data %>% rownames_to_column("cell_id") %>% as_tibble()

  write_rds(lst(expression, feature_info, sample_info), paste0(data_rds, "data.rds"))
}

if (!file.exists(paste0(data_rds, "data_filt.rds"))) {
  genesets <- read_rds("derived_files/data_genesets.rds")
  reg_entrezs <- genesets %>% filter(grepl("transcription factor", description)) %>% pull(entrezs) %>% unlist() %>% unique()

  alias2eg <- as.list(org.Mm.eg.db::org.Mm.egALIAS2EG)
  regulators <- feature_info %>% transmute(feature_id, entrez = alias2eg[feature_id]) %>% unnest(entrez) %>%
    filter(entrez %in% reg_entrezs) %>% pull(feature_id) %>% unique()

  targets_filt <- feature_info %>% arrange(desc(vst.variance.standardized)) %>% pull(feature_id) %>% head(2000)
  regulators_filt <- feature_info %>% top_n(5000, vst.variance.standardized) %>% filter(feature_id %in% regulators) %>% pull(feature_id)
  samples <- rownames(expression)

  expr_filt <- expression[, union(targets_filt, regulators_filt)]

  write_rds(lst(expr_filt, targets_filt, regulators_filt, samples, sample_info, feature_info), paste0(data_rds, "data_filt.rds"))

  rm(list = ls())
}


# SUBMIT BRED -------------------------------------------------------------
data_rds <- "derived_files/tabula_muris/data_rds/"
# list2env(read_rds(paste0(data_rds, "data.rds")), .GlobalEnv)
list2env(read_rds(paste0(data_rds, "data_filt.rds")), .GlobalEnv)


calculate_target_importance <- function(
  target_ix,
  expr,
  samples,
  regulators,
  targets,
  num_trees = 10000,
  num_variables_per_split = 100,
  num_samples_per_tree = 250,
  min_node_size = 10,
  max_depth = NULL,
  interaction_importance_filter = .01,
  sigmoid_mean = mean(expr[expr != 0]),
  sigmoid_sd = sd(expr[expr != 0])
) {

  target <- targets[[target_ix]]

  regs <- setdiff(regulators, target)

  target_expr <- scale(expr[,target])[,1]

  data <- data.frame(
    PREDICT = target_expr,
    as.matrix(expr[,regs]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  rm(expr)

  rf2 <- rangercase::ranger(
    data = data,
    dependent.variable.name = "PREDICT",
    verbose = TRUE,
    num.threads = 1,
    importance = "permutation",
    mtry = num_variables_per_split,
    num.trees = num_trees,
    min.node.size = min_node_size,
    sample.fraction = .5,
    max.depth = max_depth,
    local.importance = TRUE
  )

  cat("Selecting most important regulators\n")

  imp <- tibble(
    feature_id = names(rf2$variable.importance),
    importance = rf2$variable.importance,
    effect = ifelse(importance == 0, 0, rf2$variable.importance.cor),
    ix = seq_along(feature_id)
  ) %>%
    arrange(desc(importance))

  impf <- imp %>% filter(importance >= interaction_importance_filter)

  limp <- Matrix::t(rf2$variable.importance.casewise[impf$feature_id, , drop = FALSE])

  # downscale importance if regulator is not expressed
  # ... it can't be regulating anything if it is not expressed ...
  expr_reg_sc <- stats::pnorm(as.matrix(data[rownames(limp),colnames(limp)]), mean = sigmoid_mean, sd = sigmoid_sd)
  limp_sc <- limp * expr_reg_sc

  limp_df <-
    left_join(
      limp %>%
        reshape2::melt(varnames = c("cell_id", "regulator"), value.name = "importance"),
      limp_sc %>%
        reshape2::melt(varnames = c("cell_id", "regulator"), value.name = "importance_sc"),
      by = c("cell_id", "regulator")
    ) %>%
    as_tibble()

  lst(
    importance = imp %>% transmute(
      regulator = factor(feature_id, levels = regulators),
      target = factor(target, levels = targets),
      importance,
      effect
    ),
    importance_sc = limp_df %>% transmute(
      cell_id = factor(as.character(cell_id), levels = samples),
      regulator = factor(as.character(regulator), levels = regulators),
      target = factor(target, levels = targets),
      importance,
      importance_sc
    )
  )
}

num_trees = 10000
num_variables_per_split = 100
num_samples_per_tree = 1000
max_depth = 10
min_node_size = 10
interaction_importance_filter = .01
sigmoid_mean = mean(expr_filt@x)
sigmoid_sd = sd(expr_filt@x)

qsub_handle <- qsub::qsub_lapply(
  X = seq_along(targets_filt),
  qsub_config = qsub::override_qsub_config(
    max_wall_time = NULL,
    memory = "10G",
    name = "bred",
    wait = FALSE,
    stop_on_error = FALSE,
    remove_tmp_folder = FALSE,

  ),
  qsub_packages = c("dynutils", "dplyr", "purrr", "magrittr", "tibble"),
  qsub_environment = c("x"),
  FUN = calculate_target_importance,
  # pass data and other parameters
  expr = expr_filt,
  samples = samples,
  regulators = regulators_filt,
  targets = targets_filt,
  num_trees = num_trees,
  num_variables_per_split = num_variables_per_split,
  num_samples_per_tree = num_samples_per_tree,
  min_node_size = min_node_size,
  max_depth = max_depth,
  interaction_importance_filter = interaction_importance_filter,
  sigmoid_mean = sigmoid_mean,
  sigmoid_sd = sigmoid_sd
)
dest_dir <- "derived_files/tabula_muris/"
write_rds(qsub_handle, paste0(dest_dir, "qsub_handle.rds"))

qsub_handle <- read_rds(paste0(dest_dir, "qsub_handle.rds"))

grn <- qsub::qsub_retrieve(
  qsub_handle,
  wait = "just_do_it"
)
write_rds(grn, paste0(dest_dir, "grn.rds"))

grn <- read_rds(paste0(dest_dir, "grn.rds"))



# remove unfinished executions
grn[map_lgl(grn, ~ length(.) == 1 && is.na(.))] <- NULL

# return combined results
importance <- grn %>% map_df("importance") %>% filter(importance > .1) %>% mutate(i = row_number(), name = paste0(regulator, "->", target))
write_tsv(importance, "~/importance.tsv")
importance_sc <- grn %>%
  map_df("importance_sc") %>%
  filter(importance > .1) %>%
  inner_join(importance %>% select(-importance, -name), by = c("regulator", "target"))

imp_sc_mat <- Matrix::sparseMatrix(
  i = importance_sc$cell_id %>% as.integer,
  j = importance_sc$i,
  x = importance_sc$importance,
  dims = c(length(samples), nrow(importance)),
  dimnames = list(levels(importance_sc$cell_id), importance$name)
)

dimred <- dyndimred::dimred_landmark_mds(imp_sc_mat)
qplot(sample_info$tissue_tSNE_1, sample_info$tissue_tSNE_2)
qplot(dimred[,1], dimred[,2], col = sample_info$tissue)
g1 <- dynplot::plot_dimred(model, grouping = dataset$grouping) + ggtitle("orig dimred")
g2 <- dynplot::plot_dimred(model, dimred = dimred, grouping = dataset$grouping) + ggtitle("new dimred")
patchwork::wrap_plots(g1, g2, nrow = 1)
dynplot::plot_heatmap(model, imp_sc_mat, features_oi = 100, grouping = dataset$grouping)

dimred2 <- dyndimred::dimred_landmark_mds(Matrix::t(imp_sc_mat), ndim = 10)
cl <- kmeans(dimred2, centers = 8)

rgl::plot3d(dimred2, col = cl$cluster, size = 10)

plot_mat <- imp_sc_mat %>% as.matrix %>% dynutils::scale_quantile() %>% t()
lin <- dynplot:::linearise_cells(model)$progressions %>% arrange(edge_id, cumpercentage)
gaps_col <- which(diff(as.integer(lin$edge_id)) != 0)
pheatmap::pheatmap(
  plot_mat[order(cl$cluster),lin$cell_id],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  gaps_row = which(diff(sort(cl$cluster)) != 0),
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = dataset$grouping %>% enframe("name", "annoation") %>% as.data.frame() %>% column_to_rownames("name")
)

grouping <- model$progressions %>%
  mutate(group = case_when(
    from == "3" & percentage < .33 ~ "MEF",
    from == "3" & percentage < .66 ~ "early-induced",
    from == "3" ~ "induced",
    percentage < .33 ~ "induced",
    to == "4" & percentage < .66 ~ "early-neuron",
    to == "5" & percentage < .66 ~ "early-myocyte",
    to == "4" ~ "neuron",
    to == "5" ~ "myocyte"
  ))

imp_gr <-
  importance_sc %>%
  left_join(grouping %>% select(cell_id, group), by = "cell_id") %>%
  group_by(group, regulator, target) %>%
  summarise_at(c("importance", "effect"), mean) %>%
  ungroup()


imp_gr %>%
  filter(importance > 1) %>%
  write_tsv("~/interactions.tsv")

imp_gr %>%
  group_by(group) %>%
  top_n(50, importance) %>%
  write_tsv("~/interactions_top50.tsv")

plot_dimred(model, grouping = grouping %>% select(cell_id, group) %>% deframe())

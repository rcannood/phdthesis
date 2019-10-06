library(tidyverse)

dest_dir <- "derived_files/tabula_muris/"
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

cellontf <- paste0(dest_dir, "cl.obo")
if (!file.exists(cellontf)) download.file("https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl.obo", cellontf)

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
dest_dir <- "derived_files/tabula_muris/"
data_rds <- paste0(dest_dir, "data_rds/")
# list2env(read_rds(paste0(data_rds, "data.rds")), .GlobalEnv)
list2env(read_rds(paste0(data_rds, "data_filt.rds")), .GlobalEnv)

if (!file.exists(paste0(dest_dir, "qsub_handle.rds"))) {
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
}

dest_dir <- "derived_files/tabula_muris/"
qsub_handle <- read_rds(paste0(dest_dir, "qsub_handle.rds"))

grn <- qsub::qsub_retrieve(
  qsub_handle,
  wait = "just_do_it",
  post_fun = function(i, li) {
    li$importance <- li$importance %>% filter(importance > .01)
    li$importance_sc <- li$importance_sc %>% filter(importance_sc > .01)
    li
  }
)
write_rds(grn, paste0(dest_dir, "grn.rds"))

grn <- read_rds(paste0(dest_dir, "grn.rds"))

# remove unfinished executions
grn[map_lgl(grn, ~ length(.) == 1 && is.na(.))] <- NULL

# return combined results
importance <-
  grn %>%
  map_df("importance") %>%
  filter(importance > .1) %>%
  mutate(
    i = row_number(),
    name = paste0(regulator, "->", target)
  )
ggplot(importance) +
  geom_point(aes(effect, importance))
write_tsv(importance, "~/importance.tsv")
importance_sc <- grn %>%
  map_df("importance_sc") %>%
  filter(importance_sc > .01) %>%
  inner_join(importance %>% select(-importance, -name), by = c("regulator", "target"))

rm(grn)
gc()

list2env(read_rds(paste0(data_rds, "data_filt.rds")), .GlobalEnv)
imp_sc_mat <- Matrix::sparseMatrix(
  i = importance_sc$cell_id %>% as.integer,
  j = importance_sc$i,
  x = importance_sc$importance_sc,
  dims = c(length(samples), nrow(importance)),
  dimnames = list(levels(importance_sc$cell_id), importance$name)
)


# CELL ONTOLOGY PREP
library(ontologyIndex)
cell_ont <-
  get_ontology(cellontf) %>%
  with(tibble(id, name, parents, children, ancestors, obsolete)) %>%
  filter(!obsolete, grepl("CL:", id))

gr <- cell_ont %>% select(from = id, to = children) %>% unnest(to) %>%
  igraph::graph_from_data_frame()

dis <- igraph::distances(gr, v = "CL:0000000")
depth <- dis[1,] %>% enframe("id", "depth")

co_upstream <-
  cell_ont %>%
  select(cell_ontology_id = id, upstream = ancestors) %>%
  unnest(upstream) %>%
  filter(cell_ontology_id != upstream) %>%
  left_join(depth %>% select(upstream = id, depth), by = c("upstream"))

# clustering
dimred <- dyndimred::dimred_landmark_mds(imp_sc_mat, ndim = 10, distance_method = "spearman")
dimred_umap <- dyndimred::dimred_umap(dimred, pca_components = NULL, n_neighbors = 30)

df <- data.frame(
  # dimred,
  dimred_umap,
  sample_info
)
ggplot(df, aes(comp_1, comp_2, col = tissue)) + geom_point()

k <- 100

km <- stats::kmeans(dimred, centers = k)
clus <- km$cluster

samdf <-
  sample_info %>%
  transmute(cell_id, cell_ontology_id, cluster = clus) %>%
  filter(!is.na(cell_ontology_id))

TOT <- nrow(samdf)

test <-
  samdf %>%
  group_by(cluster) %>%
  mutate(CLUSPOS = n()) %>%
  ungroup() %>%
  left_join(co_upstream %>% rename(ontology = upstream, upstream_depth = depth), by = "cell_ontology_id") %>%
  left_join(depth %>% select(cell_ontology_id = id, orig_depth = depth), by = c("cell_ontology_id")) %>%
  mutate(diff = abs(orig_depth - upstream_depth)) %>%
  group_by(ontology) %>%
  mutate(ONTPOS = n()) %>%
  group_by(cluster, ontology) %>%
  summarise(
    CLUSPOS = CLUSPOS[[1]],
    ONTPOS = ONTPOS[[1]],
    TP = n(),
    mean_diff_depth = mean(diff)
  ) %>%
  ungroup() %>%
  mutate(
    FN = ONTPOS - TP,
    FP = CLUSPOS - TP,
    TN = TOT - TP - FN - FP,
    TPR = TP / ONTPOS,
    TNR = TN / (TN + FP),
    PPV = TP / (TP + FP),
    NPV = TN / (TN + FN),
    FNR = FN / (FN + TP),
    FPR = FP / (FP + TN),
    FDR = FP / (FP + TP),
    FOR = FN / (FN + TN),
    F12 = dynutils::calculate_harmonic_mean(cbind(PPV, TPR), weights = c(2, 1)),
    F1 = 2 * TP / (2 * TP + FP + FN),
    MMC = sqrt(PPV * TPR * TNR * NPV) - sqrt(FDR * FNR * FPR * FOR)
  ) %>%
  left_join(cell_ont %>% select(ontology = id, name), by = "ontology")
ggplot(test) + geom_point(aes(mean_diff_depth, F12))

labels <-
  test %>%
  arrange(desc(F12)) %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup()

cell_names <- unique(labels$name) %>% sort()
col_names <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(length(cell_names))

labels <- labels %>%
  mutate(col = col_names[match(name, cell_names)])

df <- data.frame(
  dimred_umap,
  sample_info
) %>%
  mutate(
    cluster = clus
  ) %>%
  left_join(labels %>% select(cluster, name), by = "cluster")
labeldf <- df %>% group_by(name) %>% summarise_at(c("comp_1", "comp_2"), mean)
g <- ggplot(df, aes(comp_1, comp_2, col = name)) +
  geom_point(size = .5) +
  # geom_text(aes(label = name), labeldf, colour = "black", size = 4, fontface = "bold") +
  geom_label(aes(label = name), labeldf, size = 4) +
  theme_bw() +
  coord_equal() +
  scale_colour_manual(values = setNames(col_names, cell_names))
g
ggsave(paste0(dest_dir, "plot_umap.pdf"), g, width = 15, height = 13)

impsc2 <-
  importance_sc %>%
  mutate(cluster = as.vector(kmfit[cell_id])) %>%
  left_join(labels %>% select(cluster, ontology, name, col), by = "cluster") %>%
  group_by(name, i, regulator, target) %>%
  summarise(
    importance = mean(importance),
    importance_sc = mean(importance_sc),
    effect = effect[[1]],
    col = col[[1]]
  ) %>%
  ungroup() %>%
  select(name, regulator, target, everything())
# impsc2f <- impsc2 %>% filter(importance_sc > 5)
impsc2 %>% filter(importance_sc > 1) %>% group_by(name) %>% summarise(n = n()) %>% as.data.frame() %>% arrange(desc(n))
impsc2f <- impsc2 %>% group_by(name) %>% arrange(desc(importance_sc)) %>% slice(1:50) %>% ungroup() %>% filter(importance_sc > .5)
impsc2f
write_tsv(impsc2f, "~/aaa_impsc2.tsv")

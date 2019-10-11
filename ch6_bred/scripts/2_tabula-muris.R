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
    # num_samples_per_tree = 250,
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
      # mtry = num_variables_per_split,
      num.trees = num_trees,
      min.node.size = min_node_size,
      # sample.fraction = .5,
      # max.depth = max_depth,
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

  num_trees = 1000
  num_variables_per_split = 100
  max_depth = 20
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
      execute_before = "#$ -l h=!prismcls05"
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
    li$importance_sc <- li$importance_sc %>% filter(importance_sc > .01) %>%
      inner_join(li$importance %>% select(regulator, target), by = c("regulator", "target"))
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
  mutate(
    i = row_number(),
    name = paste0(regulator, "->", target)
  )
ggplot(importance) +
  geom_point(aes(effect, importance))
write_tsv(importance, "~/importance.tsv")
importance_sc <- grn %>%
  map_df("importance_sc") %>%
  inner_join(importance %>% select(-importance, -name), by = c("regulator", "target"))

rm(grn)
gc()

list2env(read_rds(paste0(data_rds, "data_filt.rds")), .GlobalEnv)
list2env(read_rds("derived_files/cell_ontology.rds"), .GlobalEnv)

# clustering
imp_sc_mat <- Matrix::sparseMatrix(
  i = importance_sc$cell_id %>% as.integer,
  j = importance_sc$i,
  x = importance_sc$importance_sc,
  dims = c(length(samples), nrow(importance)),
  dimnames = list(samples, importance$name)
)
dimred <- dyndimred::dimred_landmark_mds(imp_sc_mat, ndim = 20, distance_method = "spearman")
rm(imp_sc_mat)

knn <- RANN::nn2(dimred, k = 100)
knn$nn.dists <- knn$nn.dists[,-1]
knn$nn.idx <- knn$nn.idx[,-1]

knndf <-
  inner_join(
    reshape2::melt(knn$nn.idx, varnames = c("i", "nn"), value.name = "j"),
    reshape2::melt(knn$nn.dists, varnames = c("i", "nn"), value.name = "dist"),
    by = c("i", "nn")
  ) %>%
  as_tibble() %>%
  mutate(i2 = pmin(i, j), j2 = pmax(i, j)) %>%
  select(i = i2, j = j2, dist) %>%
  unique() %>%
  arrange(dist) %>%
  mutate(weight = 1 / dist)

gr <- igraph::graph_from_data_frame(knndf %>% select(i, j, dist), vertices = seq_len(nrow(dimred)), directed = FALSE)
cl <- igraph::cluster_louvain(gr)
clus <- cl$membership
dimred2 <- igraph::layout_with_fr(gr)
rownames(dimred2) <- rownames(dimred)
colnames(dimred2) <- paste0("comp_", seq_len(ncol(dimred2)))


knns <- pbapply::pblapply(
  seq_len(12),
  cl = 1,
  function(repi) {
    cix <- sample.int(ncol(knn$nn.idx), floor(.9 * ncol(knn$nn.idx)))
    knndf <-
      inner_join(
        reshape2::melt(knn$nn.idx[,cix], varnames = c("i", "nn"), value.name = "j"),
        reshape2::melt(knn$nn.dists[,cix], varnames = c("i", "nn"), value.name = "dist"),
        by = c("i", "nn")
      ) %>%
      as_tibble() %>%
      mutate(i2 = pmin(i, j), j2 = pmax(i, j)) %>%
      select(i = i2, j = j2, dist) %>%
      unique() %>%
      arrange(dist) %>%
      mutate(weight = 1 / dist)

    gr <- igraph::graph_from_data_frame(knndf %>% select(i, j, dist), vertices = seq_len(nrow(dimred)), directed = FALSE)
    dimred2 <- igraph::layout_with_fr(gr)
    rownames(dimred2) <- rownames(dimred)
    colnames(dimred2) <- paste0("comp_", seq_len(ncol(dimred2)))

    lst(dimred2)
  }
)

dimred_fr <- map(knns, "dimred2") %>% do.call(cbind, .) %>% dynutils::scale_uniform()

dimred_fr2 <- dyndimred::dimred_pca(dimred_fr)
plot(dimred_fr2, col = clus)

dimred_umap <- dyndimred::dimred_umap(dimred_fr, pca_components = NULL, n_neighbors = 50)
plot(dimred_umap, col = clus)

df <- data.frame(
  dimred_umap,
  # dimred2,
  # dimred_fr2,
  sample_info,
  clus
)



`%s/%` <- function(x, y) ifelse(y == 0, 1, x / y)
ks <- seq(10, 30)
backgr <- pbapply::pblapply(ks, cl = 1, function(k) {
  out <- map_df(seq_len(100), function(i) {
    km <- stats::kmeans(dimred_fr, centers = k)
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
      left_join(cell_ont_up %>% rename(ontology = upstream), by = "cell_ontology_id") %>%
      group_by(ontology) %>%
      mutate(ONTPOS = n()) %>%
      group_by(cluster, ontology) %>%
      summarise(
        CLUSPOS = CLUSPOS[[1]],
        ONTPOS = ONTPOS[[1]],
        TP = n()
      ) %>%
      ungroup() %>%
      mutate(
        FN = ONTPOS - TP,
        FP = CLUSPOS - TP,
        TN = TOT - TP - FN - FP,
        TPR = TP %s/% (TP + FN),
        TNR = TN %s/% (TN + FP),
        PPV = TP %s/% (TP + FP),
        NPV = TN %s/% (TN + FN),
        FNR = FN %s/% (FN + TP),
        FPR = FP %s/% (FP + TN),
        FDR = FP %s/% (FP + TP),
        FOR = FN %s/% (FN + TN),
        F1 = (2 * TP) %s/% (2 * TP + FP + FN)
      )
  })

  out %>% select(-cluster:-TN) %>% map(ecdf)
})

outs <- pbapply::pblapply(rep(seq_along(ks), 10), cl = 1, function(ki) {
  k <- ks[[ki]]

  km <- stats::kmeans(dimred_fr, centers = k)
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
    left_join(cell_ont_up %>% rename(ontology = upstream), by = "cell_ontology_id") %>%
    group_by(ontology) %>%
    mutate(ONTPOS = n()) %>%
    group_by(cluster, ontology) %>%
    summarise(
      CLUSPOS = CLUSPOS[[1]],
      ONTPOS = ONTPOS[[1]],
      TP = n()
    ) %>%
    ungroup() %>%
    left_join(cell_ont %>% select(ontology = id, name), by = "ontology") %>%
    mutate(
      FN = ONTPOS - TP,
      FP = CLUSPOS - TP,
      TN = TOT - TP - FN - FP,
      TPR = TP %s/% (TP + FN),
      TNR = TN %s/% (TN + FP),
      PPV = TP %s/% (TP + FP),
      NPV = TN %s/% (TN + FN),
      FNR = FN %s/% (FN + TP),
      FPR = FP %s/% (FP + TN),
      FDR = FP %s/% (FP + TP),
      FOR = FN %s/% (FN + TN),
      F1 = (2 * TP) %s/% (2 * TP + FP + FN)
    )

  ecdfs <- backgr[[ki]]
  for (ne in names(ecdfs)) {
    test[[paste0("p_", ne)]] <- ecdfs[[ne]](test[[ne]])
  }

  labels <-
    test %>%
    arrange(desc(F1)) %>%
    group_by(cluster) %>%
    slice(1) %>%
    ungroup()

  summ <- labels %>% summarise_if(is.numeric, dynutils::calculate_geometric_mean) %>% mutate(k, ki)

  lst(clus, test, labels, summ, km)
})

summ <- map_df(outs, "summ") %>% mutate(row = row_number())
ggplot(summ) + geom_point(aes(k, F1))
ggplot(summ) + geom_point(aes(k, p_F1))

out <- outs[[summ %>% arrange(desc(p_F1), desc(F1)) %>% pull(row) %>% first()]]
out$summ
km <- out$km
labels <- out$labels
clus <- out$clus


labels <- labels %>%
  group_by(name) %>%
  mutate(name2 = if (n() == 1) name else paste0(name, " #", row_number())) %>%
  ungroup() %>%
  select(-name) %>%
  rename(name = name2)


cell_names <- unique(labels$name) %>% sort()
col_names <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(length(cell_names))

labels <- labels %>%
  mutate(col = col_names[match(name, cell_names)])

df <- data.frame(
  dimred_fr2,
  sample_info
) %>%
  mutate(
    cluster = clus
  ) %>%
  left_join(labels %>% select(cluster, name), by = "cluster")


# labeldf <- df %>% group_by(cluster) %>% summarise(comp_1 = mean(comp_1), comp_2 = mean(comp_2), name = name[[1]])
labeldf <- df %>% group_by(name) %>% summarise(comp_1 = mean(comp_1), comp_2 = mean(comp_2))
g <- ggplot(df, aes(comp_1, comp_2, col = name)) +
  geom_point(size = .5) +
  shadowtext::geom_shadowtext(aes(label = name), labeldf, bg.colour = "white", size = 5) +
  theme_bw() +
  coord_equal() +
  scale_colour_manual(values = setNames(col_names, cell_names))
g
ggsave(paste0(dest_dir, "plot_umap.pdf"), g, width = 15, height = 13)

gc()

impsc2 <-
  importance_sc %>%
  mutate(cluster = as.vector(clus[cell_id])) %>%
  left_join(labels %>% transmute(cluster, name = factor(name, levels = cell_names)), by = "cluster") %>%
  group_by(name, i, regulator, target) %>%
  summarise(
    importance = mean(importance),
    importance_sc = mean(importance_sc),
    effect = effect[[1]]
  ) %>%
  ungroup() %>%
  mutate(col = col_names[name]) %>%
  select(regulator, name, target, everything())
# impsc2f <- impsc2 %>% filter(importance_sc > 5)
impsc2f <- impsc2 %>% group_by(name) %>% arrange(desc(importance_sc)) %>% slice(1:50) %>% ungroup() %>% filter(importance_sc > 1)
impsc2f %>% group_by(name) %>% summarise(n = n()) %>% arrange(desc(n))
write_tsv(impsc2f, paste0(dest_dir, "aaa_impsc2.tsv"))











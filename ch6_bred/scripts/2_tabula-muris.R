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

if (!file.exists(paste0(dest_dir, "grn.rds"))) {
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
}
grn <- read_rds(paste0(dest_dir, "grn.rds"))




# remove unfinished executions
grn[map_lgl(grn, ~ length(.) == 1 && is.na(.))] <- NULL

# return combined results
importance <-
  grn %>%
  map_df("importance") %>%
  arrange(desc(importance)) %>%
  mutate(
    interaction_id = forcats::fct_inorder(paste0(regulator, "->", target))
  )
importance_sc <-
  grn %>%
  map_df("importance_sc") %>%
  rename(sample_id = cell_id) %>%
  left_join(importance %>% select(regulator, target, effect, interaction_id), by = c("regulator", "target"))

rm(grn)
gc()



if (!file.exists(paste0(dest_dir, "dimred_and_clustering.rds"))) {
  dimred_out <- bred:::dimred_and_cluster(importance_sc, knn = 100, use_scaled_imp = TRUE)

  plot(dimred_out$dimred_fr, col = dimred_out$cluster)
  # apply gng
  # gng_fit <- gng::gng(dimred_out$dimred_lmds, max_nodes = 20)
  # dimred_out$dimred_fr <- gng_fit$space_proj
  # dimred_out$cluster <- match(gng_fit$clustering, gng_fit$nodes$name) %>% setNames(names(gng_fit$clustering))

  cluster <- dimred_out$cluster

  # annotate the clusters
  cellont <- read_rds("derived_files/cell_ontology.rds")
  sample_groupings_cellont <-
    sample_info %>%
    select(sample_id = cell_id, co_id = cell_ontology_id) %>%
    left_join(
      cellont$cell_ont %>% select(co_id = id, co_anc = ancestors) %>% filter(co_id %in% unique(sample_info$cell_ontology_id)) %>% unnest(co_anc),
      by = "co_id"
    ) %>%
    select(-co_id) %>%
    left_join(cellont$cell_ont %>% select(co_anc = id, group_id = name), by = "co_anc") %>%
    transmute(sample_id, group_type = "cell_ontology", group_id)

  sample_groupings_clin <-
    sample_info %>%
    select(sample_id = cell_id, mouse.id, mouse.sex) %>%
    gather(group_type, group_id, -sample_id)

  sample_groupings <- bind_rows(
    sample_groupings_cellont,
    sample_groupings_clin
  )
  relabel_out <- bred::label_clusters(sample_groupings, cluster, arrange_fun = function(df) df %>% arrange(desc(F1)))

  relabel_out$labels <- relabel_out$labels %>% mutate(name = paste0(name, ifelse(PPV < .5, " *", "")))

  relabel_out$labels %>% print(n = 100)

  # look at cluster annotation
  group_type_palette <-
    unique(sample_groupings$group_type) %>%
    {setNames(RColorBrewer::brewer.pal(length(.), "Dark2"), .)}
  pdf(paste0(dest_dir, "cluster_labelling.pdf"), width = 12, height = 6)
  tryCatch({
    for (cli in relabel_out$labels$cluster) {
      test1c <- relabel_out$test %>%
        filter(cluster == cli) %>%
        mutate(x = forcats::fct_reorder(paste0(group_type, ": ", group_id), F1)) %>%
        arrange(x) %>%
        top_n(30, F1) %>%
        gather(metric, value, PPV, TPR, F1)
      g <- ggplot(test1c, aes(x, value, fill = group_type)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        coord_flip() +
        labs(x = NULL, title = paste0("Cluster ", cli, ", ", relabel_out$labels$name[[cli]])) +
        theme(legend.position = "bottom") +
        scale_fill_manual(values = group_type_palette) +
        scale_x_discrete(labels = test1c$group_id) +
        facet_wrap(~metric, nrow = 1)

      print(g)
    }
  }, finally = {
    dev.off()
  })

  # save dimred and clustering
  write_rds(c(dimred_out, relabel_out, lst(sample_groupings)), paste0(dest_dir, "dimred_and_clustering.rds"))
}
# list2env(c(dimred_out, relabel_out), .GlobalEnv)
list2env(read_rds(paste0(dest_dir, "dimred_and_clustering.rds")), .GlobalEnv)



palette <- labels %>% select(name, col) %>% deframe()

sample_info_df <- data.frame(
  dimred_fr[sample_info$cell_id,],
  sample_info,
  cluster = cluster[sample_info$cell_id],
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  left_join(labels, by = "cluster") %>%
  as_tibble()
labeldf <- sample_info_df %>% group_by(name) %>% summarise(comp_1 = mean(comp_1), comp_2 = mean(comp_2))

# make heatmap
tab <- sample_info_df %>%
  group_by(cell_ontology_class, name) %>%
  summarise(n = n()) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  reshape2::acast(name~cell_ontology_class, value.var = "pct", fill = 0)
tab <- tab[, order(apply(tab, 2, which.max))]
gaps_col <- which(diff(apply(tab, 2, which.max)) != 0)
pheatmap::pheatmap(tab, gaps_col = gaps_col, angle_col = 315, cluster_rows = FALSE, cluster_cols = FALSE, filename = paste0(dest_dir, "cluster_heatmap.pdf"), width = 16, height = 8)

# save visualisation
g <-
  ggplot(sample_info_df, aes(comp_1, comp_2, col = name)) +
  geom_point(size = .5) +
  shadowtext::geom_shadowtext(aes(label = name), labeldf, bg.colour = "white", size = 5) +
  dynplot::theme_graph() +
  coord_equal() +
  scale_colour_manual(values = palette) +
  theme(legend.position = "none")
g
ggsave(paste0(dest_dir, "plot_fr.pdf"), g, width = 10, height = 10)



g <-
  ggplot() +
  # geom_point(size = .5) +
  geom_segment(aes(x = 1, xend = 2, y = name, yend = name, col = name), labeldf, size = 2) +
  theme_bw() +
  scale_colour_manual(values = palette) +
  labs(colour = "Group") +
  theme(legend.position = "bottom") +
  guides(
    colour = guide_legend(nrow = 6)
  )
g
ggsave(paste0(dest_dir, "legend.pdf"), g, width = 15, height = 6)

clusn <- table(cluster)
imp_grouped <-
  importance_sc %>%
  mutate(cluster = as.vector(cluster[sample_id])) %>%
  left_join(labels %>% select(cluster, name), by = "cluster") %>%
  group_by(name, interaction_id, regulator, target) %>%
  summarise(
    importance = sum(importance) / clusn[[cluster[[1]]]],
    importance_sc = sum(importance_sc) / clusn[[cluster[[1]]]],
    effect = effect[[1]]
  ) %>%
  ungroup() %>%
  left_join(labels %>% select(name, col), by = "name") %>%
  select(source = regulator, name, target, everything()) %>%
  arrange(desc(importance))

imp_grouped_f <- imp_grouped %>% group_by(name) %>% slice(1:50) %>% ungroup()
imp_grouped_f %>% group_by(name) %>% summarise(n = n()) %>% arrange(desc(n)) %>% print(n = Inf)
write_tsv(imp_grouped_f, paste0(dest_dir, "grouped_interactions.tsv"))

write_tsv(imp_grouped %>% filter(grepl("B cell", name)) %>% group_by(name) %>% slice(1:50) %>% ungroup(), paste0(dest_dir, "grouped_interactions_b.tsv"))
write_tsv(imp_grouped %>% filter(grepl("T cell", name)) %>% group_by(name) %>% slice(1:100) %>% ungroup(), paste0(dest_dir, "grouped_interactions_t.tsv"))

gng_fit <- gng::gng(dimred_lmds)

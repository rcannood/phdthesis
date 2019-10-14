library(tidyverse)

tcga_folder <- "~/Workspace/tcga/"
data_dir <- "derived_files/tcga/"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

tcga_data_file <- paste0(tcga_folder, "data.rds")
data_file <- paste0(data_dir, "data.rds")

if (!file.exists(tcga_data_file)) {
  sample_info <- read_tsv(paste0(tcga_folder, "sample_sheet.tsv"))
  colnames(sample_info) <- colnames(sample_info) %>% tolower() %>% gsub(" ", "_", .)

  # preproc metadata
  sample_info$project_id <- gsub(",.*", "", sample_info$project_id)
  sample_info$sample_type <- gsub(",.*", "", sample_info$sample_type)
  sample_info$sample_id <- gsub(",.*", "", sample_info$sample_id)
  sample_info <- sample_info %>%
    filter(!sample_type %in% c("Control Analyte", "Additional Metastatic", "FFPE Scrolls", "Slides")) %>%
    mutate(
      sample_group = case_when(
        sample_type == "Blood Derived Normal" ~ "Normal Blood",
        sample_type == "Solid Tissue Normal" ~ "Normal Tissue",
        grepl("Tumor", sample_type) ~ "Cancer Tumor",
        sample_type == "Additional - New Primary" ~ "Cancer Tumor",
        grepl("Metastatic", sample_type) ~ "Cancer Metastatic",
        grepl("Cancer.*Bone Marrow", sample_type) ~ "Cancer BM",
        grepl("Cancer.*Blood", sample_type) ~ "Cancer Blood"
      ),
      file = paste0(tcga_folder, "data/", file_id, "/", file_name)
    ) %>%
    group_by(project_id, sample_group) %>%
    mutate(
      id = paste0(project_id, "_", tolower(gsub(" ", "-", sample_group)), "_sample", seq_len(n()))
    ) %>%
    ungroup() %>%
    select(id, project_id, sample_group, everything())

  # check classes
  table(sample_info$project_id) %>% sort
  table(sample_info$sample_type) %>% sort
  table(sample_info$sample_group) %>% sort
  table(sample_info$sample_type, sample_info$sample_group)

  # check genes and rename genes^
  feature_info <-
    read_tsv(sample_info$file[[1]], col_names = c("id", "values"), col_types = cols(.default = "d", id = "c")) %>%
    transmute(orig_id = id, ensembl_gene_id = gsub("\\..*", "", id))

  ensembl <- biomaRt::useMart(
    "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "http://mar2015.archive.ensembl.org"
  )
  gene_symbols <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene"),
    filters = "ensembl_gene_id",
    values = feature_info$ensembl_gene_id,
    mart = ensembl
  ) %>%
    as_tibble()

  gene_symbols_gath <-
    gene_symbols %>%
    gather(type, value, hgnc_symbol, entrezgene) %>%
    filter(!is.na(value), value != "") %>%
    group_by(ensembl_gene_id, type) %>%
    summarise(value = list(value)) %>%
    spread(type, value)

  feature_info <-
    left_join(
      feature_info,
      gene_symbols_gath,
      by = "ensembl_gene_id"
    ) %>%
    mutate(
      symbol = ifelse(map_lgl(hgnc_symbol, is.null), ensembl_gene_id, map_chr(hgnc_symbol, ~ ifelse(length(.) > 0, first(.), NA_character_)))
    ) %>%
    group_by(symbol) %>%
    mutate(
      id = if (n() > 1) paste0(symbol, "_", row_number()) else symbol
    ) %>%
    ungroup() %>%
    select(id, ensembl_gene_id, symbol, everything())

  # load in all files
  tmp_dir <- "~/tcga_tmp/"
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  by <- round(nrow(sample_info) / 100)
  start_ix <- seq(1, nrow(sample_info), by = by)
  end_ix <- pmin(nrow(sample_info), start_ix + by - 1)
  mat_files <- paste0(tmp_dir, "mat_", seq_along(start_ix), ".rds")
  exists <- file.exists(mat_files)

  assertthat::assert_that(
    !any(duplicated(sample_info$id)),
    !any(duplicated(feature_info$id))
  )

  rm <- pbapply::pblapply(which(!exists), cl = 8, function(i) {
    mat_file <- mat_files[[i]]
    if (!file.exists(mat_file)) {
      cat("Progress ", i, "/", length(start_ix), "\n", sep = "")
      six <- start_ix[[i]]
      eix <- end_ix[[i]]
      ixs <- seq(six, eix)
      mats <- map(ixs, function(j) {
        mat <-
          read_tsv(
            sample_info$file[[j]],
            col_names = c("feature_id", "counts"),
            col_types = cols(counts = "d", feature_id = "c"),
            progress = FALSE
          ) %>%
          column_to_rownames("feature_id") %>%
          as.matrix()
        if (!all(rownames(mat) == feature_info$orig_id)) {
          if (!all(rownames(mat) %in% feature_info$orig_id)) {
            stop("Rowname mismatch in sample ", i)
          }
          mat <- mat[feature_info$orig_id, ]
        }
        rownames(mat) <- NULL
        mat
      }) %>% do.call(cbind, .)

      rownames(mats) <- feature_info$id
      colnames(mats) <- sample_info$id[ixs]

      write_rds(mats, mat_file)

      NULL
    }
  })

  counts <- map(seq_along(start_ix), function(i) {
    mat_file <- paste0(tmp_dir, "mat_", i, ".rds")
    read_rds(mat_file)
  }) %>%
    do.call(cbind, .) %>%
    t()

  # if (file.exists(paste0(tcga_folder, "rnaseq.rds"))) {
  #   unlink(tmp_dir, recursive = TRUE)
  # }

  gc()

  write_rds(lst(counts, sample_info, feature_info, gene_symbols), tcga_data_file)
}

if (!file.exists(data_file)) {
  list2env(read_rds(tcga_data_file), .GlobalEnv)
  expr <- log2(counts + 1)

  rm(counts)

  means <- colMeans(expr)
  vars <- apply(expr, 2, var)
  minexpr <- colMeans(expr > 2)

  expr <- expr[, minexpr >= .01]

  feature_info <- feature_info %>% slice(match(colnames(expr), id))

  genesets <- read_rds("derived_files/data_genesets_human.rds")
  reg_entrezs <- genesets %>% filter(grepl("transcription factor", description)) %>% pull(entrezs) %>% unlist() %>% unique()
  regulators <- feature_info %>% select(id, entrezgene) %>% unnest(entrezgene) %>% filter(entrezgene %in% reg_entrezs) %>% pull(id) %>% unique()

  vars <- apply(expr, 2, var)
  targets <- names(sort(vars, decreasing = TRUE))

  samples <- sample_info$id

  assertthat::assert_that(
    !any(duplicated(sample_info$id)),
    !any(duplicated(feature_info$id)),
    !any(duplicated(rownames(expr))),
    !any(duplicated(colnames(expr)))
  )

  write_rds(lst(expr, sample_info, feature_info, regulators, targets, samples), data_file)
}





# SUBMIT BRED -------------------------------------------------------------
if (!file.exists(paste0(data_dir, "qsub_handle.rds"))) {
  list2env(read_rds(data_file), .GlobalEnv)

  calculate_target_importance <- function(
    target_ix,
    expr,
    samples,
    regulators,
    targets,
    num_trees = 10000,
    num_variables_per_split = 100,
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
  min_node_size = 20
  interaction_importance_filter = .01
  sigmoid_mean = mean(expr)
  sigmoid_sd = sd(expr)

  regulators_filt <- regulators
  targets_filt <- targets[1:2000]
  expr_filt <- expr[, union(targets_filt, regulators_filt)]
  samples_filt <- rownames(expr_filt)

  x <- 1
  qsub_handle <- qsub::qsub_lapply(
    X = seq_along(targets_filt),
    qsub_config = qsub::override_qsub_config(
      max_wall_time = "12:00:00",
      memory = "10G",
      name = "bred",
      wait = FALSE,
      stop_on_error = FALSE,
      remove_tmp_folder = FALSE,
      execute_before = "#$ -l h=!prismcls08",
    ),
    qsub_packages = c("dynutils", "dplyr", "purrr", "magrittr", "tibble"),
    qsub_environment = c("x"),
    FUN = calculate_target_importance,
    # pass data and other parameters
    expr = expr_filt,
    samples = samples_filt,
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
  write_rds(qsub_handle, paste0(data_dir, "qsub_handle.rds"))

  targets_filt <- targets[2001:4000]
  expr_filt <- expr[, union(targets_filt, regulators_filt)]

  x <- 1
  qsub_handle2 <- qsub::qsub_lapply(
    X = seq_along(targets_filt),
    qsub_config = qsub::override_qsub_config(
      max_wall_time = "12:00:00",
      memory = "10G",
      name = "bred",
      wait = FALSE,
      stop_on_error = FALSE,
      remove_tmp_folder = FALSE#,
      # execute_before = "#$ -l h=!prismcls08",
    ),
    qsub_packages = c("dynutils", "dplyr", "purrr", "magrittr", "tibble"),
    qsub_environment = c("x"),
    FUN = calculate_target_importance,
    # pass data and other parameters
    expr = expr_filt,
    samples = samples_filt,
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
  write_rds(qsub_handle2, paste0(data_dir, "qsub_handle2.rds"))
}


if (!file.exists(paste0(data_dir, "grn.rds"))) {
  qsub_handle1 <- read_rds(paste0(data_dir, "qsub_handle.rds"))

  grn1 <- qsub::qsub_retrieve(
    qsub_handle1,
    wait = "just_do_it",
    post_fun = function(i, li) {
      li$importance <- li$importance %>% filter(importance > .01)
      li$importance_sc <- li$importance_sc %>% filter(importance_sc > .01) %>%
        inner_join(li$importance %>% select(regulator, target), by = c("regulator", "target"))
      li
    }
  )

  qsub_handle2 <- read_rds(paste0(data_dir, "qsub_handle2.rds"))

  grn2 <- qsub::qsub_retrieve(
    qsub_handle2,
    wait = "just_do_it",
    post_fun = function(i, li) {
      li$importance <- li$importance %>% filter(importance > .01)
      li$importance_sc <- li$importance_sc %>% filter(importance_sc > .01) %>%
        inner_join(li$importance %>% select(regulator, target), by = c("regulator", "target"))
      li
    }
  )

  grn <- c(grn1, grn2)
  rm(grn1, grn2, qsub_handle1, qsub_handle2)

  write_rds(grn, paste0(data_dir, "grn.rds"))
}

list2env(read_rds(data_file), .GlobalEnv)

if (!file.exists(paste0(data_dir, "importances.rds"))) {
  grn <- read_rds(paste0(data_dir, "grn.rds"))

  # remove unfinished executions
  grn[map_lgl(grn, ~ length(.) == 1 && is.na(.))] <- NULL

  # return combined results
  importance <-
    grn %>%
    map_df(function(li) {
      li$importance %>% mutate(
        regulator = factor(regulator, levels = colnames(expr)),
        target = factor(target, levels = colnames(expr))
      )
    }) %>%
    arrange(desc(importance)) %>%
    mutate(
      i = row_number(),
      name = forcats::fct_inorder(paste0(regulator, "->", target)),
    )
  ggplot(importance) +
    geom_point(aes(effect, importance))

  write_tsv(importance %>% slice(1:850), paste0(data_dir, "interactions.tsv"))

  importance_sc <-
    grn %>%
    map_df(function(li) {
      li$importance_sc %>% mutate(
        regulator = factor(regulator, levels = colnames(expr)),
        target = factor(target, levels = colnames(expr))
      )
    }) %>%
    inner_join(importance %>% select(-importance, -name), by = c("regulator", "target")) %>%
    arrange(desc(importance))

  rm(grn)
  gc()

  write_rds(lst(importance, importance_sc), paste0(data_dir, "importances.rds"))
}


list2env(read_rds(paste0(data_dir, "importances.rds")), .GlobalEnv)

if (!file.exists(paste0(data_dir, "dimred_and_clustering.rds"))) {
  # calculate LMDS dimred
  samples <- levels(importance_sc$cell_id)
  imp_sc_mat <- Matrix::sparseMatrix(
    i = importance_sc$cell_id %>% as.integer,
    j = importance_sc$i,
    x = importance_sc$importance,
    dims = c(length(samples), nrow(importance)),
    dimnames = list(samples, importance$name)
  )

  dimred <- dyndimred::dimred_landmark_mds(imp_sc_mat, ndim = 20, distance_method = "spearman")
  rm(imp_sc_mat)
  gc()

  # compute knn
  knn <- RANN::nn2(dimred, k = 101)
  knn$nn.dists <- knn$nn.dists[,-1]
  knn$nn.idx <- knn$nn.idx[,-1]

  # perform louvain clustering and FR dimred
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
  dimred_fr <- igraph::layout_with_fr(gr)


  return_rotation_mat <- function(theta) {
    theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  rot_mat <- return_rotation_mat(-90)
  colnames(rot_mat) <- colnames(dimred_fr)
  plot(dimred_fr %*% rot_mat)

  dimred_fr <- dimred_fr %*% rot_mat

  rownames(dimred_fr) <- rownames(dimred)
  colnames(dimred_fr) <- paste0("comp_", seq_len(ncol(dimred_fr)))

  # annotate the clusters
  proj_annot <-
    read_tsv(paste0(data_dir, "projects.tsv")) %>%
    mutate(
      wheretype = ifelse(!is.na(ontology) & !is.na(type), paste0(gsub(",.*", "", where), " ", type), NA_character_)
    ) %>%
    select(-ontology, -name_long) %>%
    rename(project_id = project) %>%
    gather(group, value, -project_id) %>%
    na.omit() %>%
    mutate(value = strsplit(value, ", ")) %>%
    unnest(value)

  samdf <-
    sample_info %>%
    transmute(id, project_id, cluster = clus)

  TOT <- nrow(samdf)

  `%s/%` <- function(x, y) ifelse(y == 0, 1, x / y)
  test <-
    samdf %>%
    group_by(cluster) %>%
    mutate(CLUSPOS = n()) %>%
    ungroup() %>%
    left_join(proj_annot, by = "project_id") %>%
    group_by(group, value) %>%
    mutate(GRPOS = n()) %>%
    group_by(cluster, group, value) %>%
    summarise(
      CLUSPOS = CLUSPOS[[1]],
      GRPOS = GRPOS[[1]],
      TP = n()
    ) %>%
    ungroup() %>%
    mutate(
      FN = GRPOS - TP,
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
      F1 = (2 * PPV * TPR) %s/% (PPV + TPR)
    )

  labels <-
    test %>%
    arrange(desc(PPV * TPR)) %>%
    group_by(cluster) %>%
    slice(1) %>%
    ungroup() %>%
    rename(name = value) %>%
    group_by(name) %>%
    mutate(name2 = factor(if (n() == 1) name else paste0(name, " #", row_number()))) %>%
    ungroup() %>%
    select(-name) %>%
    rename(name = name2) %>%
    arrange(name) %>%
    mutate(col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(n())) %>%
    arrange(cluster)

  palette <- labels %>% select(name, col) %>% deframe()

  # look at cluster annotation
  clus <- cl$membership
  tab <- as.matrix(table(
    paste0(clus, " ", labels$name[clus]),
    sample_info$project_id
  ))
  tab <- sweep(tab, 1, rowSums(tab), "/")
  tab <- tab[order(apply(tab, 1, max), decreasing = TRUE), ]
  tab <- tab[, order(apply(tab, 2, which.max))]
  pheatmap::pheatmap(tab, angle_col = 315, cluster_rows = FALSE, cluster_cols = FALSE, filename = paste0(data_dir, "cluster_heatmap.pdf"), width = 16, height = 8)

  # save dimred and clustering
  write_rds(lst(knn, dimred_fr, clus, labels, palette), paste0(data_dir, "dimred_and_clustering.rds"))
}

list2env(read_rds(paste0(data_dir, "dimred_and_clustering.rds")), .GlobalEnv)

# save visualisation
df <- data.frame(
  dimred_fr,
  sample_info,
  cluster = clus
) %>%
  left_join(labels, by = "cluster")
labeldf <- df %>% group_by(name) %>% summarise(comp_1 = mean(comp_1), comp_2 = mean(comp_2))

g <-
  ggplot(df, aes(comp_1, comp_2, col = name)) +
  geom_point(size = .5) +
  shadowtext::geom_shadowtext(aes(label = name), labeldf, bg.colour = "white", size = 5) +
  theme_bw() +
  coord_equal() +
  scale_colour_manual(values = palette)
g
ggsave(paste0(data_dir, "plot_fr.pdf"), g, width = 15, height = 13)

g <-
  ggplot() +
  # geom_point(size = .5) +
  geom_segment(aes(x = 1, xend = 2, y = name, yend = name, col = name), labeldf, size = 2) +
  theme_bw() +
  scale_colour_manual(values = palette) +
  labs(colour = "Group")
g
ggsave(paste0(data_dir, "legend.pdf"), g, width = 6, height = 6)


imp_grouped <-
  importance_sc %>%
  mutate(cluster = as.vector(clus[cell_id])) %>%
  left_join(labels %>% select(cluster, name), by = "cluster") %>%
  group_by(name, i, regulator, target) %>%
  summarise(
    importance = mean(importance),
    importance_sc = mean(importance_sc),
    effect = effect[[1]]
  ) %>%
  ungroup() %>%
  left_join(labels %>% select(name, col), by = "name") %>%
  select(regulator, name, target, everything()) %>%
  arrange(desc(importance))

imp_grouped_f <- imp_grouped %>% group_by(name) %>% slice(1:50) %>% ungroup() %>% filter(importance > .15)
imp_grouped_f %>% group_by(name) %>% summarise(n = n()) %>% arrange(desc(n)) %>% print(n = Inf)
write_tsv(imp_grouped_f, paste0(data_dir, "grouped_interactions.tsv"))


# grpdm <- Matrix::sparseMatrix(
#   i = as.integer(imp_grouped$name),
#   j = imp_grouped$i,
#   x = imp_grouped$importance,
#   dims = c(nrow(labels), nrow(importance)),
#   dimnames = list(
#     labels$name,
#     importance$name
#   )
# )
# sim <- dynutils::calculate_similarity(grpdm, method = "spearman")
# diag(sim) <- NA
# pheatmap::pheatmap(sim, cutree_rows = 4, cutree_cols = 4)

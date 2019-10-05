library(tidyverse)
library(dyno)

if (!file.exists("derived_files/dahlin2.rds")) {
  # download and preprocess dataset
  dest_dir <- "derived_files/GSE107727/"
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  # COUNTS ------------------------------------------------------------------

  # file <- paste0(dest_dir, "data.tar")
  # url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE107727&format=file"
  # if (!file.exists(file)) {
  #   download.file(url, file)
  # }
  # if (!file.exists(paste0(dest_dir, "GSM2877127_SIGAB1_counts.txt.gz"))) {
  #   untar(file, exdir = dest_dir)
  # }
  #
  # count_files <- list.files(dest_dir, pattern = "counts\\.txt\\.gz", full.names = TRUE)
  # counts <-
  #   pbapply::pblapply(count_files, function(file) {
  #     gsm <- gsub(".*(GSM\\d*)_.*", "\\1", file)
  #     mat <-
  #       read.csv(file, header = TRUE, sep = "\t") %>%
  #       as.matrix %>%
  #       Matrix::Matrix(sparse = TRUE) %>%
  #       Matrix::t()
  #     rownames(mat) <- paste0(gsm, "_", rownames(mat))
  #     mat
  #   }) %>% do.call(rbind, .)
  # write_rds(counts, paste0(dest_dir, "counts.rds"), compress = "gz")
  counts <- read_rds(paste0(dest_dir, "counts.rds"))

  # GENE INFO ---------------------------------------------------------------
  ensembl <- biomaRt::useMart(
    "ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host = "http://mar2016.archive.ensembl.org"
  )
  gene_symbols <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene"),
    filters = "ensembl_gene_id",
    values = colnames(counts),
    mart = ensembl
  ) %>%
    as_tibble()

  feature_info <-
    gene_symbols %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      symbol = mgi_symbol[[1]],
      entrez = list(na.omit(entrezgene)),
      all_symbols = list(na.omit(mgi_symbol))
    ) %>%
    rename(feature_id = ensembl_gene_id)

  # CELL INFO ---------------------------------------------------------------
  geo <- GEOquery::getGEO("GSE107727", destdir = dest_dir)
  phenodata <- geo[[1]] %>%
    Biobase::phenoData() %>%
    as("data.frame")
  split <- phenodata %>% map(~length(unique(.)) == 1) %>% unlist()

  experiment_info <- phenodata[,split]
  sample_info <- phenodata[,!split]
  cell_info <-
    tibble(cell_id = rownames(counts)) %>%
    mutate(
      sample_id = gsub("_.*", "", cell_id),
      barcode = gsub(".*_", "", cell_id)
    ) %>%
    left_join(
      sample_info %>%
        select(
          sample_id = geo_accession,
          cell_type = `cell type:ch1`,
          genotype = `genotype:ch1`,
          strain = `strain:ch1`
        ),
      by = "sample_id"
    )

  write_rds(lst(counts, cell_info, feature_info, gene_symbols), "derived_files/dahlin2.rds", compress = "gz")
}

list2env(read_rds("derived_files/dahlin2.rds"), .GlobalEnv)

# cell_info <- cell_info %>% filter(genotype == "WT")
# counts <- counts[cell_info$cell_id, ]

seu <-
  Seurat::CreateSeuratObject(Matrix::t(counts), meta.data = cell_info %>% column_to_rownames("cell_id"), min.cells = 3, min.features = 200) %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData()
rm(counts)

genesets <- read_rds("derived_files/data_genesets.rds")
reg_entrezs <- genesets %>% filter(grepl("transcription factor", description)) %>% pull(entrezs) %>% unlist() %>% unique()
regulators <- gene_symbols %>% filter(entrezgene %in% reg_entrezs) %>% pull(ensembl_gene_id) %>% unique()


targets <- seu@assays$RNA@var.features
regulators <- intersect(regulators, targets)
expr <- seu@assays$RNA@data[targets, ] %>% Matrix::t() %>% as.matrix
samples <- rownames(expr)


# RUN BRED ----------------------------------------------------------------

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
  interaction_importance_filter = .01,
  sigmoid_mean = mean(expr[expr != 0]),
  sigmoid_sd = sd(expr[expr != 0]),
  num_perms = 10
) {

  target <- targets[[target_ix]]

  regs <- setdiff(regulators, target)

  target_expr <- scale(expr[,target])[,1]

  # data <- Matrix::Matrix(expr[,union(target, regs)], sparse = TRUE)
  # data[,1] <- target_expr
  # colnames(data)[[1]] <- "PREDICT"

  data <- data.frame(
    PREDICT = target_expr,
    expr[,regs],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # st1 <- system.time({
  rf2 <- ranger::ranger(
    data = data,
    dependent.variable.name = "PREDICT",
    verbose = TRUE,
    num.threads = 1,
    importance = "impurity",
    mtry = num_variables_per_split,
    num.trees = num_trees,
    min.node.size = min_node_size,
    sample.fraction = num_samples_per_tree / nrow(expr),
    keep.inbag = TRUE,
    save.memory = TRUE
  )
  # })

  # st2 <- system.time({
  # rf <- randomForest::randomForest(
  #   expr[,regs],
  #   target_expr,
  #   mtry = num_variables_per_split,
  #   ntree = num_trees,
  #   sampsize = num_samples_per_tree,
  #   nodesize = min_node_size,
  #   importance = FALSE,
  #   localImp = FALSE,
  #   keep.inbag = TRUE
  # )
  # })

  imp <-
    rf2$variable.importance %>%
    enframe("feature_id", "importance") %>%
    mutate(ix = row_number()) %>%
    # rf$importance %>%
    # as.data.frame() %>%
    # tibble::rownames_to_column("feature_id") %>%
    # tibble::as_tibble() %>%
    # mutate(ix = row_number()) %>%
    # rename(importance = `%IncMSE`) %>%
    arrange(desc(importance))

  impf <- imp %>% filter(importance >= interaction_importance_filter * 100) # if IncNodePurity
  impf$effect <- NA

  reg_check <- impf$feature_id
  limp <- matrix(0, nrow = length(samples), ncol = length(reg_check), dimnames = list(samples, reg_check))

  # inbag <- rf$inbag
  inbag <- do.call(cbind, rf2$inbag.counts)

  shuffles <- map(seq_len(num_perms), ~ sample.int(nrow(expr)))
  pred <- predict(rf2, expr, predict.all = TRUE, num.threads = 1)
  # pred_ind <- pred$individual
  pred_ind <- pred$predictions
  pred_ind[inbag > 0] <- NA
  pred_ind_sw <- sweep(pred_ind, 1, target_expr, "-")^2
  pred_agg <- rowMeans(pred_ind, na.rm = TRUE)

  expr2 <- expr
  eff_x <- eff_y <- rep(NA, nrow(expr) * num_perms)

  for (j in seq_along(reg_check)) {
    regj <- reg_check[[j]]
    cat("Running permutations for regulator ", j, " / ", length(reg_check), ": ", regj, "\n", sep = "")

    expr_regj <- expr[,regj]

    for (i in seq_len(num_perms)) {
      ix <- shuffles[[i]]

      expr2[,regj] <- expr_regj[ix]

      pred2 <- predict(rf2, expr2, predict.all = TRUE, num.threads = 1)
      # pred2_ind <- pred2$individual
      pred2_ind <- pred2$predictions
      pred2_ind[inbag > 0] <- NA

      inc_se <-
        sweep(pred2_ind, 1, target_expr, "-")^2 - pred_ind_sw

      # calculate importance
      imp <- rowMeans(inc_se, na.rm = TRUE)

      limp[, j] <- limp[, j] + imp / num_perms

      # determine effect
      effix <- seq_len(nrow(expr)) + nrow(expr) * (i - 1)
      eff_x[effix] <- expr2[,regj] - expr[,regj]
      eff_y[effix] <- rowMeans(pred2_ind, na.rm = TRUE) - pred_agg
    }
    impf$effect[[j]] <- cor(eff_x, eff_y)

    # reset expr
    expr2[,regj] <- expr_regj
  }

  # we've now computed the permutation importance
  impf$importance <- colMeans(limp)

  # downscale importance if regulator is not expressed
  # ... it can't be regulating anything if it is not expressed ...
  expr_reg_sc <- stats::pnorm(expr[,reg_check], mean = sigmoid_mean, sd = sigmoid_sd)
  limp_sc <- limp * expr_reg_sc

  limp_df <- limp_sc %>%
    reshape2::melt(varnames = c("cell_id", "regulator"), value.name = "importance") %>%
    as_tibble()

  lst(
    importance = impf %>% transmute(
      regulator = factor(feature_id, levels = regulators),
      target = factor(target, levels = targets),
      importance,
      effect
    ),
    importance_sc = limp_df %>% transmute(
      cell_id = factor(as.character(cell_id), levels = samples),
      regulator = factor(as.character(regulator), levels = regulators),
      target = factor(target, levels = targets),
      importance = importance
    )
  )
}

num_trees = 10000
num_variables_per_split = 100
num_samples_per_tree = 250
min_node_size = 10
interaction_importance_filter = .01
sigmoid_mean = mean(expr@x)
sigmoid_sd = sd(expr@x)
num_perms = 10

x <- 1
qsub_handle <- qsub::qsub_lapply(
  X = seq_along(targets),
  qsub_config = qsub::override_qsub_config(
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
  expr = expr,
  samples = samples,
  regulators = regulators,
  targets = targets,
  num_trees = num_trees,
  num_variables_per_split = num_variables_per_split,
  num_samples_per_tree = num_samples_per_tree,
  min_node_size = min_node_size,
  interaction_importance_filter = interaction_importance_filter,
  sigmoid_mean = sigmoid_mean,
  sigmoid_sd = sigmoid_sd,
  num_perms = num_perms
)
write_rds(qsub_handle, "derived_files/qsub_handle.rds")

qsub_handle <- read_rds("derived_files/qsub_handle.rds")

grn <- qsub::qsub_retrieve(
  qsub_handle,
  wait = "just_do_it"
)
write_rds(grn, "derived_files/grn.rds")






# FAILED LABELLING EXPERIMENT ---------------------------------------------


# # regulators <- dataset$feature_info %>% filter(entrez %in% reg_entrezs, feature_id %in% targets) %>% pull(feature_id)
#
#
# seu <- Seurat::RunPCA(seu, verbose = FALSE, npcs = 20)
# seu <- Seurat::RunUMAP(seu, dims = 1:20, umap.method = "uwot", n.neighbors = 7)
# seu <- Seurat::RunUMAP(seu, dims = 1:20, umap.method = "umap-learn", n.neighbors = 7)
#
# celltype_markers <-
#   tribble(
#     ~cell_type, ~markers, ~ensembl,
#     "HSC", "Procr", "ENSMUSG00000027611",
#     "Erythroid", "Gata1", "ENSMUSG00000031162",
#     "Erythroid", "Klf1", "ENSMUSG00000054191",
#     "Erythroid", "Epor", "ENSMUSG00000006235",
#     "Erythroid", "Gypa", "ENSMUSG00000051839",
#     "Erythroid", "Hba-a2", "ENSMUSG00000069917",
#     "Erythroid", "Hba-a1", "ENSMUSG00000069919",
#     "Erythroid", "Spi1", "ENSMUSG00000002111",
#     "Neutrophil", "Elane", "ENSMUSG00000020125",
#     "Neutrophil", "Cebpe", "ENSMUSG00000052435",
#     "Neutrophil", "Ctsg", "ENSMUSG00000040314",
#     "Neutrophil", "Mpo", "ENSMUSG00000009350",
#     "Neutrophil", "Gfi1", "ENSMUSG00000029275",
#     "Monocyte", "Irf8", "ENSMUSG00000041515",
#     "Monocyte", "Csf1r", "ENSMUSG00000024621",
#     "Monocyte", "Ctsg", "ENSMUSG00000040314",
#     "Monocyte", "Mpo", "ENSMUSG00000009350",
#     "Megakaryocyte", "Itga2b", "ENSMUSG00000034664",
#     "Megakaryocyte", "Pbx1", "ENSMUSG00000052534",
#     "Megakaryocyte", "Sdpr", "ENSMUSG00000045954",
#     "Megakaryocyte", "Vwf", "ENSMUSG00000001930",
#     "Basophil", "Mcpt8", "ENSMUSG00000022157",
#     "Basophil", "Prss34", "ENSMUSG00000056399",
#     "B cell", "Cd19", "ENSMUSG00000030724",
#     "B cell", "Vpreb2", "ENSMUSG00000059280",
#     "B cell", "Cd79a", "ENSMUSG00000003379",
#     "Mast cell", "Cma1", "ENSMUSG00000022225",
#     "Mast cell", "Gzmb", "ENSMUSG00000015437",
#     "Mast cell", "Kit", "ENSMUSG00000005672", # CD117/c-Kit
#     "Mast cell & Basophil", "Ms4a2", "ENSMUSG00000024680",
#     "Mast cell & Basophil", "Fcer1a", "ENSMUSG00000005339",
#     "Mast cell & Basophil", "Cpa3", "ENSMUSG00000001865",
#     "Mast cell & Basophil", "Enpp3", "ENSMUSG00000019989", # CD203c human
#   )
#
# Seurat::DimPlot(seu, reduction = "umap",pt.size = 0.5, label = TRUE, repel = TRUE)
#
#
# dimred <- seu@reductions$umap@cell.embeddings
# colnames(dimred) <- paste0("comp_", seq_len(ncol(dimred)))
# plot(dimred)
#
# df <- data.frame(
#   seu@meta.data %>% rownames_to_column("cell_id"),
#   dimred,
#   dynutils::scale_quantile(t(as.matrix(seu@assays$RNA@data[celltype_markers$ensembl,])))
#   # dynutils::scale_quantile(log2(as.matrix(counts[rownames(dimred),celltype_markers$ensembl])+1))
# ) %>%
#   sample_n(n()) %>%
#   as_tibble()
#
# pt_size <- .5
# patchwork::wrap_plots(
#   ggplot(df) + geom_point(aes(comp_1, comp_2, colour = cell_type), size = pt_size) + theme_bw() + scale_colour_brewer(palette = "Set1") + coord_equal(),
#   ggplot(df) + geom_point(aes(comp_1, comp_2, colour = genotype), size = pt_size) + theme_bw() + scale_colour_brewer(palette = "Set2") + coord_equal(),
#   ggplot(df) + geom_point(aes(comp_1, comp_2, colour = strain), size = pt_size) + theme_bw() + scale_colour_brewer(palette = "Set3") + coord_equal(),
#   ggplot(df) + geom_point(aes(comp_1, comp_2, colour = sample_id), size = pt_size) + theme_bw() + scale_colour_brewer(palette = "Dark2") + coord_equal(),
#   nrow = 2
# )
#
# df_expr <- df %>%
#   gather(ensembl, value, starts_with("ENS")) %>%
#   inner_join(celltype_markers %>% rename(cell_type2 = cell_type), by = "ensembl") %>%
#   group_by(cell_id, cell_type2) %>%
#   summarise(value = mean(value), comp_1 = comp_1[[1]], comp_2 = comp_2[[1]]) %>%
#   ungroup()
#
# ggplot(df_expr) +
#   geom_point(aes(comp_1, comp_2, colour = value)) +
#   facet_wrap(~cell_type2) +
#   theme_bw() +
#   viridis::scale_color_viridis() +
#   coord_equal()

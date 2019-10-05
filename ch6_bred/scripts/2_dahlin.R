library(tidyverse)
library(dyno)

if (!file.exists("derived_files/dahlin.rds")) {
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

  num_reads_per_gene <- Matrix::colSums(counts)
  qplot(log10(num_reads_per_gene + 1)) + geom_vline(xintercept = log10(10 + 1))
  num_reads_per_cell <- Matrix::rowSums(counts)
  qplot(log10(num_reads_per_cell + 1)) + geom_vline(xintercept = log10(2000 + 1))
  num_genes_expressed <- Matrix::rowSums(counts > 0)
  qplot(log10(num_genes_expressed + 1)) + geom_vline(xintercept = log10(2000 + 1))

  counts <- counts[
    num_reads_per_cell > 2000 & num_genes_expressed > 2000,
    num_reads_per_gene > 10
  ]


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

  # SAVE DATASET ------------------------------------------------------------
  expression <- counts / Matrix::rowMeans(counts) * 1000000
  expression@x <- log2(expression@x + 1)
  dataset <-
    wrap_expression(
      counts = counts,
      expression = expression,
      cell_info = cell_info,
      feature_info = feature_info,
      feature_mapper = gene_symbols
    )
  rm(list = setdiff(ls(), "dataset"))
  write_rds(dataset, "derived_files/dahlin.rds", compress = "gz")
}
dataset <- read_rds("derived_files/dahlin.rds")

celltype_markers <-
  tribble(
    ~cell_type, ~markers, ~ensembl,
    "HSC", "Procr", "ENSMUSG00000027611",
    "Erythroid", "Gata1", "ENSMUSG00000031162",
    "Erythroid", "Klf1", "ENSMUSG00000054191",
    "Erythroid", "Epor", "ENSMUSG00000006235",
    "Erythroid", "Gypa", "ENSMUSG00000051839",
    "Erythroid", "Hba-a2", "ENSMUSG00000069917",
    "Erythroid", "Hba-a1", "ENSMUSG00000069919",
    "Erythroid", "Spi1", "ENSMUSG00000002111",
    "Neutrophil", "Elane", "ENSMUSG00000020125",
    "Neutrophil", "Cebpe", "ENSMUSG00000052435",
    "Neutrophil", "Ctsg", "ENSMUSG00000040314",
    "Neutrophil", "Mpo", "ENSMUSG00000009350",
    "Neutrophil", "Gfi1", "ENSMUSG00000029275",
    "Monocyte", "Irf8", "ENSMUSG00000041515",
    "Monocyte", "Csf1r", "ENSMUSG00000024621",
    "Monocyte", "Ctsg", "ENSMUSG00000040314",
    "Monocyte", "Mpo", "ENSMUSG00000009350",
    "Megakaryocyte", "Itga2b", "ENSMUSG00000034664",
    "Megakaryocyte", "Pbx1", "ENSMUSG00000052534",
    "Megakaryocyte", "Sdpr", "ENSMUSG00000045954",
    "Megakaryocyte", "Vwf", "ENSMUSG00000001930",
    "Basophil", "Mcpt8", "ENSMUSG00000022157",
    "Basophil", "Prss34", "ENSMUSG00000056399",
    "B cell", "Cd19", "ENSMUSG00000030724",
    "B cell", "Vpreb2", "ENSMUSG00000059280",
    "B cell", "Cd79a", "ENSMUSG00000003379",
    "Mast cell", "Cma1", "ENSMUSG00000022225",
    "Mast cell", "Gzmb", "ENSMUSG00000015437",
    "Mast cell", "Kit", "ENSMUSG00000005672", # CD117/c-Kit
    "Mast cell & Basophil", "Ms4a2", "ENSMUSG00000024680",
    "Mast cell & Basophil", "Fcer1a", "ENSMUSG00000005339",
    "Mast cell & Basophil", "Cpa3", "ENSMUSG00000001865",
    "Mast cell & Basophil", "Enpp3", "ENSMUSG00000019989", # CD203c human
  )
dataset$counts[,celltype_markers$ensembl]

cell_info <- dataset$cell_info %>% filter(genotype == "WT")
counts <- dataset$counts[cell_info$cell_id, ]
rm(dataset)

gc()
seu <-
  Seurat::CreateSeuratObject(Matrix::t(counts), meta.data = cell_info %>% column_to_rownames("cell_id"), min.cells = 3, min.features = 200) %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData()
seu <- Seurat::RunPCA(seu, verbose = FALSE, npcs = 20)
seu <- Seurat::RunUMAP(seu, dims = 1:20, umap.method = "uwot", n.neighbors = 7)
Seurat::DimPlot(seu, reduction = "umap",pt.size = 0.5, label = TRUE, repel = TRUE)

# pca.results <- irlba::irlba(A = dataset$expression, nv = 20)
# total.variance <- sum(Seurat:::SparseRowVar(Matrix::t(dataset$expression), display_progress = F))
# feature.loadings <- pca.results$v
# sdev <- pca.results$d/sqrt(max(1, ncol(dataset$expression) - 1))
# predr <- pca.results$u %*% diag(pca.results$d)
# # predr <- lmds::lmds(dataset$expression, distance_method = "spearman", num_landmarks = 1000, ndim = 20)
# # predr <- dyndimred::dimred_pca(dataset$expression, ndim = 20)
# dimred <- dyndimred::dimred_umap(predr, ndim = 3, pca_components = NULL, n_neighbors = 7)
# # cell.embeddings %>% as.data.frame() %>% gather(col, value) %>% ggplot() + geom_violin(aes(col, value))
# # dimred <- dyndimred::dimred_umap(dataset$expression, distance_method = "cosine")
# # dimred <- lmds::lmds(dataset$expression, distance_method = "spearman", num_landmarks = 1000, ndim = 3)

dimred <- seu@reductions$umap@cell.embeddings
colnames(dimred) <- paste0("comp_", seq_len(ncol(dimred)))
plot(dimred)

df <- data.frame(
  cell_info,
  dimred,
  dynutils::scale_quantile(log2(as.matrix(counts[,celltype_markers$ensembl])+1))
) %>%
  sample_n(10000) %>%
  as_tibble()

patchwork::wrap_plots(
  ggplot(df) + geom_point(aes(comp_1, comp_2, colour = cell_type)) + theme_bw() + scale_colour_brewer(palette = "Set1") + coord_equal(),
  ggplot(df) + geom_point(aes(comp_1, comp_2, colour = genotype)) + theme_bw() + scale_colour_brewer(palette = "Set2") + coord_equal(),
  ggplot(df) + geom_point(aes(comp_1, comp_2, colour = strain)) + theme_bw() + scale_colour_brewer(palette = "Set3") + coord_equal(),
  ggplot(df) + geom_point(aes(comp_1, comp_2, colour = sample_id)) + theme_bw() + scale_colour_brewer(palette = "Dark2") + coord_equal(),
  nrow = 2
)

df_expr <- df %>%
  gather(ensembl, value, starts_with("ENS")) %>%
  inner_join(celltype_markers %>% rename(cell_type2 = cell_type), by = "ensembl") %>%
  group_by(cell_id, cell_type2) %>%
  summarise(value = mean(value), comp_1 = comp_1[[1]], comp_2 = comp_2[[1]]) %>%
  ungroup()

ggplot(df_expr) +
  geom_point(aes(comp_1, comp_2, colour = value)) +
  facet_wrap(~cell_type2) +
  theme_bw() +
  viridis::scale_color_viridis() +
  coord_equal()

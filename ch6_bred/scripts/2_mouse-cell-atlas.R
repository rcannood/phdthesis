library(tidyverse)
library(rlang)

# download and preprocess dataset
dest_dir <- "derived_files/mouse-cell-atlas/"
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

# cell info
cell_assignment_file <- paste0(dest_dir, "cellassignment.csv")
if (!file.exists(cell_assignment_file)) download.file("https://ndownloader.figshare.com/files/11083451?private_link=865e694ad06d5857db4b", cell_assignment_file)
all_cell_info <-
  read_csv(cell_assignment_file) %>%
  rename(cell_id = Cell.name, cluster_id = ClusterID, tissue = Tissue, batch = Batch, barcode = Cell.Barcode, group_id = Annotation) %>%
  mutate(
    cell_type = group_id %>% str_extract("^[^\\(_]*"),
    subtype = group_id %>% str_extract("_([^\\()]*)\\(") %>% str_replace_all("[_\\(]", "")
  )
all_cell_info %>% select(cell_type) %>% unique %>% arrange(cell_type) %>% write_tsv(paste0(dest_dir, "/clean_cell_types.tsv"))
all_cell_info %>% select(cell_type, tissue) %>% unique %>% arrange(cell_type) %>% write_tsv(paste0(dest_dir, "/clean_cell_types2.tsv"))
# counts
counts_dir <- paste0(dest_dir, "counts/")
counts_files <- list.files(counts_dir, pattern = ".txt.gz", full.names = TRUE, recursive = TRUE)

if (length(counts_files) == 0) {
  counts_zip_file <- paste0(dest_dir, "count_batchremove.zip")
  download.file("https://ndownloader.figshare.com/files/10756795?private_link=865e694ad06d5857db4b", counts_zip_file)

  dir.create(counts_dir, showWarnings = FALSE, recursive = FALSE)
  unzip(counts_zip_file, exdir = counts_dir)
  counts_files <- list.files(counts_dir, pattern = ".txt.gz", full.names = TRUE, recursive = TRUE)
}

counts <- pbapply::pblapply(counts_files, function(file) {
  rds <- gsub("\\.txt\\.gz$", ".rds", file)
  if (!file.exists(rds)) {
    mat <- read.csv(file, sep = " ") %>%
      as.matrix() %>%
      Matrix::Matrix(sparse = TRUE) %>%
      Matrix::t()
    write_rds(mat, rds, compress = "gz")
  } else {
    read_rds(rds)
  }
})
colnamess <- map(counts, colnames) %>% unlist() %>% unique() %>% sort()

for (i in seq_along(counts)) {
  cat(i, "/", length(counts), "\n", sep = "")
  co <- counts[[i]]
  em <- Matrix::Matrix(0, sparse = TRUE, ncol = 1, nrow = nrow(co), dimnames = list(rownames(co), "empty"))
  coc <- cbind(co, em)

  ix <- match(colnamess, colnames(co)) %|% ncol(coc)
  newc <- coc[, ix]
  colnames(newc) <- colnamess
  counts[[i]] <- newc
}

counts <- do.call(rbind, counts)

write_rds(counts, paste0(dest_dir, "counts.rds"))

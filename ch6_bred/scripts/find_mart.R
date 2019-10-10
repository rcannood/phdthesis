

# MOUSE -------------------------------------------------------------------
archs <- biomaRt::listEnsemblArchives()
gene_symbols <- pbapply::pblapply(seq_len(nrow(archs)), function(i) {
  tryCatch({
    ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host=archs$url[[i]])
    attrs <- biomaRt::listAttributes(ensembl)
    entrez_attr_name <- ifelse("entrezgene" %in% attrs$name, "entrezgene", "entrezgene_id")
    gene_symbols <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "mgi_symbol", entrez_attr_name),
      filters = "ensembl_gene_id",
      values = colnames(counts),
      mart = ensembl
    )
    gene_symbols %>%
      rename(entrez_gene_id = !!entrez_attr_name) %>%
      as_tibble()
  }, error = function(e) {
    NULL
  })
})


archs %>%
  mutate(
    matching_genes = map_int(gene_symbols, function(x) if(is.data.frame(x)) length(unique(x$ensembl_gene_id)) else 0L),
    num_genes = ncol(counts)
  ) %>%
  filter(num_genes == matching_genes)





# HUMAN -------------------------------------------------------------------
archs <- biomaRt::listEnsemblArchives()
gene_symbols <- pbapply::pblapply(seq_len(nrow(archs)), function(i) {
  tryCatch({
    ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host=archs$url[[i]])
    attrs <- biomaRt::listAttributes(ensembl)
    entrez_attr_name <- ifelse("entrezgene" %in% attrs$name, "entrezgene", "entrezgene_id")
    gene_symbols <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", entrez_attr_name),
      filters = "ensembl_gene_id",
      values = feature_info$feature_id %>% gsub("\\..*$", "", .),
      mart = ensembl
    )
    gene_symbols %>%
      rename(entrez_gene_id = !!entrez_attr_name) %>%
      as_tibble()
  }, error = function(e) {
    NULL
  })
})


archs %>%
  mutate(
    matching_genes = map_int(gene_symbols, function(x) if(is.data.frame(x)) length(unique(x$ensembl_gene_id)) else 0L),
    num_genes = nrow(feature_info)
  ) %>%
  filter(num_genes == matching_genes)

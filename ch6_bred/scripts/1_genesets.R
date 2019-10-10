library(tidyverse)

## CREATE GO tibble
go.entrezs <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
ontologies <- annotate::Ontology(names(go.entrezs))
go.entrezs <- go.entrezs[ontologies %in% c("BP", "MF")]
go.entrezs <- lapply(go.entrezs, function(go) na.omit(unique(go)))
go.ids <- names(go.entrezs)
go.terms <- annotate::Term(GO.db::GOTERM[go.ids])
go <- tibble(database = "GO", id = go.ids, description = go.terms, entrezs = go.entrezs)

## CREATE KEGG OBJECTS
pws <- pbapply::pblapply(names(KEGGREST::keggList("pathway", "mmu")), function(n) KEGGREST::keggGet(n)[[1]])
kegg.entrezs <- lapply(pws, function(pw) {
  genes <- pw$GENE
  if (length(genes) > 0) unique(genes[seq(1, length(genes)-1, by=2)])
  else genes
})
kegg.ids <- map_chr(pws, ~ .$ENTRY)
kegg.names <- map_chr(pws, ~ .$NAME)
kegg <- tibble(
  database = "KEGG",
  id = kegg.ids,
  description = kegg.names %>% gsub(" - [^-]*", "", .),
  entrezs = kegg.entrezs
)

## CREATE REACTOME OBJECTS
PATHID2NAME <- as.list(reactome.db::reactomePATHID2NAME)
PATHID2NAME <- PATHID2NAME[grepl("Mus musculus", PATHID2NAME)]
EXTID2PATHID <- as.list(reactome.db::reactomeEXTID2PATHID)
reactome.entrezs <- pbapply::pblapply(names(PATHID2NAME), cl = 8, function(id) names(EXTID2PATHID)[sapply(EXTID2PATHID, function(x) id %in% x)])
reactome <- tibble(
  database = "reactome",
  id = names(PATHID2NAME),
  description = PATHID2NAME %>% gsub("^Mus musculus: ", "", .),
  entrezs = reactome.entrezs
)

# genesets db
genesets <- bind_rows(
  go,
  kegg,
  reactome
)

write_rds(genesets, "derived_files/data_genesets_mouse.rds")









## CREATE GO tibble
go.entrezs <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
ontologies <- annotate::Ontology(names(go.entrezs))
go.entrezs <- go.entrezs[ontologies %in% c("BP", "MF")]
go.entrezs <- lapply(go.entrezs, function(go) na.omit(unique(go)))
go.ids <- names(go.entrezs)
go.terms <- annotate::Term(GO.db::GOTERM[go.ids])
go <- tibble(database = "GO", id = go.ids, description = go.terms, entrezs = go.entrezs)

## CREATE KEGG OBJECTS
pws <- pbapply::pblapply(names(KEGGREST::keggList("pathway", "hsa")), function(n) KEGGREST::keggGet(n)[[1]])
kegg.entrezs <- lapply(pws, function(pw) {
  genes <- pw$GENE
  if (length(genes) > 0) unique(genes[seq(1, length(genes)-1, by=2)])
  else genes
})
kegg.ids <- map_chr(pws, ~ .$ENTRY)
kegg.names <- map_chr(pws, ~ .$NAME)
kegg <- tibble(
  database = "KEGG",
  id = kegg.ids,
  description = kegg.names %>% gsub(" - [^-]*", "", .),
  entrezs = kegg.entrezs
)

## CREATE REACTOME OBJECTS
PATHID2NAME <- as.list(reactome.db::reactomePATHID2NAME)
PATHID2NAME <- PATHID2NAME[grepl("Homo sapiens", PATHID2NAME)]
EXTID2PATHID <- as.list(reactome.db::reactomeEXTID2PATHID)
reactome.entrezs <- pbapply::pblapply(names(PATHID2NAME), cl = 8, function(id) names(EXTID2PATHID)[sapply(EXTID2PATHID, function(x) id %in% x)])
reactome <- tibble(
  database = "reactome",
  id = names(PATHID2NAME),
  description = PATHID2NAME %>% gsub("^Homo sapiens: ", "", .),
  entrezs = reactome.entrezs
)

# genesets db
genesets <- bind_rows(
  go,
  kegg,
  reactome
)

write_rds(genesets, "derived_files/data_genesets_human.rds")

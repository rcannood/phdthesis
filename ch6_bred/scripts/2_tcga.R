library(tidyverse)

tcga_folder <- "~/Workspace/tcga/"
data_dir <- "derived_files/tcga/"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

tcga_data_file <- paste0(tcga_folder, "data.rds")
data_file <- paste0(data_dir, "data.rds")

if (!file.exists(tcga_data_file)) {
  sample_sheet <- read_tsv(paste0(tcga_folder, "sample_sheet.tsv"))
  colnames(sample_sheet) <- colnames(sample_sheet) %>% tolower() %>% gsub(" ", "_", .)

  # preproc metadata
  sample_sheet$project_id <- gsub(",.*", "", sample_sheet$project_id)
  sample_sheet$sample_type <- gsub(",.*", "", sample_sheet$sample_type)
  sample_sheet$sample_id <- gsub(",.*", "", sample_sheet$sample_id)
  sample_sheet <- sample_sheet %>%
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
      old_id = paste0(project_id, "_", tolower(gsub(" ", "-", sample_group)), "_sample", seq_len(n()))
    ) %>%
    ungroup() %>%
    filter(!duplicated(file_id))

  assertthat::assert_that(!any(duplicated(sample_sheet$file_id)))

  # add clinical case id
  mjs <- c(
    jsonlite::read_json("/home/rcannood/Workspace/tcga/raw/F_metadata.cart.2019-10-08.json"),
    jsonlite::read_json("/home/rcannood/Workspace/tcga/raw/M_metadata.cart.2019-10-08.json")
  )

  meta <- pbapply::pblapply(mjs, function(mji) {
    if (length(mji$associated_entities) > 0) {
      tibble(
        file_id = mji$file_id,
        file_name = mji$file_name,
        clinical.case_id = mji$associated_entities[[1]]$case_id
      )
    } else {
      NULL
    }
  }) %>% bind_rows() %>% unique()
  sample_sheet2 <- sample_sheet %>%
    left_join(meta, by = c("file_id", "file_name"))
  assertthat::assert_that(!any(duplicated(sample_sheet2$file_id)))

  # add clinical data
  cjs <- c(
    jsonlite::read_json("/home/rcannood/Workspace/tcga/raw/F_clinical/clinical.json"),
    jsonlite::read_json("/home/rcannood/Workspace/tcga/raw/M_clinical/clinical.json")
  )

  clinical <- pbapply::pblapply(seq_along(cjs), cl = 8, function(i) {
    jsi <- cjs[[i]]
    if (length(jsi$diagnoses) > 0) {
      jsi$diagnosis <- jsi$diagnoses[[1]]
      jsi$diagnosis$treatments <- NULL
    }

    if (length(jsi$exposures) > 0) {
      jsi$exposure <- jsi$exposures[[1]]
    }

    li <- list(case_id = jsi$case_id)
    for (n in intersect(c("diagnosis", "demographic", "exposure"), names(jsi))) {
      vals <- jsi[[n]]
      li[paste0(n, ".", names(vals))] <- map(vals, ~ . %||% NA)
    }
    as_tibble(li)
  }) %>% bind_rows()

  clinical2 <- clinical %>%
    transmute(
      clinical.case_id = case_id,
      gender = demographic.gender,
      race = demographic.race,
      ethnicity = demographic.ethnicity,
      vital_status = demographic.vital_status,
      last_known_disease_status = diagnosis.last_known_disease_status,
      days_to_birth = demographic.days_to_birth,
      days_to_diagnosis = diagnosis.days_to_diagnosis,
      days_to_last_follow_up = diagnosis.days_to_last_follow_up,
      primary_diagnosis = tolower(diagnosis.primary_diagnosis),
      tissue_or_organ = tolower(diagnosis.tissue_or_organ_of_origin)
    )

  sample_sheet3 <-
    sample_sheet2 %>%
    left_join(clinical2, by = "clinical.case_id")

  assertthat::assert_that(!any(duplicated(sample_sheet3$file_id)))

  if (!file.exists(paste0(tcga_folder, "/maps/map_out.tsv"))) {
    proj_meta <- read_tsv(paste0(tcga_folder, "/raw/project_meta.tsv"))
    tabn <-
      clinical %>%
      transmute(case_id, primary_diagnosis = tolower(diagnosis.primary_diagnosis), tissue_or_organ = tolower(diagnosis.tissue_or_organ_of_origin)) %>%
      left_join(meta, by = "case_id") %>%
      left_join(sample_sheet %>% select(project_id, file_id, file_name) %>% unique(), by = c("file_id", "file_name")) %>%
      group_by(primary_diagnosis, tissue_or_organ, project_id) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      left_join(proj_meta, by = "project_id") %>%
      mutate_all(~ ifelse(is.na(.), "", .))
    write_tsv(tabn, paste0(tcga_folder, "/maps/map_in.tsv"))
    # process manually <scream>
  }
  nciont <- read_tsv(paste0(tcga_folder, "/maps/map_out.tsv")) %>% select(primary_diagnosis, tissue_or_organ, nci_ont) %>% unique()

  sample_info <-
    sample_sheet3 %>%
    left_join(nciont, by = c("primary_diagnosis", "tissue_or_organ")) %>%
    select(-file_name, -data_category, -data_type, -case_id, -sample_id, -sample_type) %>%
    rename(id = file_id)

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

  rm(sample_sheet, sample_sheet2, sample_sheet3, mjs, cjs, clinical, clinical2, ensembl, gene_symbols_gath)

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

  rm <- pbapply::pblapply(which(!exists), cl = 1, function(i) {
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
list2env(read_rds(data_file), .GlobalEnv)

if (!file.exists(paste0(data_dir, "qsub_handle.rds"))) {
  regulators_filt <- regulators
  targets_filt <- targets[1:5000]
  expr_filt <- expr[, union(targets_filt, regulators_filt)]
  samples_filt <- rownames(expr_filt)

  qsub_handle <- bred::qsub_submit_bred(
    expr_filt,
    samples = samples_filt,
    regulators = regulators_filt,
    targets = targets_filt,
    num_trees = 1000,
    min_node_size = 10,
    sigmoid_mean = mean(expr),
    sigmoid_sd = sd(expr)
  )
  write_rds(qsub_handle, paste0(data_dir, "qsub_handle.rds"))
}

if (!file.exists(paste0(data_dir, "importances.rds"))) {
  qsub_handle <- read_rds(paste0(data_dir, "qsub_handle.rds"))

  out <- qsub::qsub_retrieve(
    qsub_handle,
    wait = "just_do_it",
    post_fun = function(i, li) {
      li$importance <- li$importance %>% filter(importance > .01)
      li$importance_sc <- li$importance_sc %>% filter(importance_sc > .01) %>%
        inner_join(li$importance %>% select(regulator, target), by = c("regulator", "target"))
      li
    }
  )

  write_rds(out, paste0(data_dir, "importances.rds"))
}

list2env(read_rds(paste0(data_dir, "importances.rds")), .GlobalEnv)

# fixes for old results
assertthat::assert_that(all(levels(importance_sc$cell_id) %in% sample_info$old_id))

importance <- importance %>% select(-i) %>% rename(interaction_id = name)
importance_sc <-
  importance_sc %>%
  rename(old_id = cell_id) %>%
  left_join(sample_info %>% transmute(sample_id = factor(id), old_id = factor(old_id, levels = levels(importance_sc$cell_id))), "old_id") %>%
  filter(!is.na(sample_id)) %>%
  select(-i) %>%
  left_join(importance %>% select(regulator, target, interaction_id), by = c("regulator", "target"))

if (!file.exists(paste0(data_dir, "dimred_and_clustering.rds"))) {
  dimred_out <- bred::dimred_and_cluster(importance_sc)

  # annotate the clusters
  nciont <- read_rds("derived_files/nci_ontology.rds")
  sample_groupings_nci <-
    sample_info %>%
    select(sample_id = id, nci_ont) %>%
    left_join(
      nciont %>% select(nci_ont = id, nci_anc = ancestors) %>% filter(nci_ont %in% unique(sample_info$nci_ont)) %>% unnest(nci_anc),
      by = "nci_ont"
    ) %>%
    select(-nci_ont) %>%
    left_join(nciont %>% select(nci_anc = id, group_id = name), by = "nci_anc") %>%
    transmute(sample_id, group_type = "nci_ontology", group_id)

  sample_groupings_clin <-
    sample_info %>%
    left_join(read_tsv(paste0(data_dir, "projects.tsv")) %>% select(project_id = project, project_tag = name_short), by = "project_id") %>%
    select(sample_id = id, project_tag, sample_group, gender, ethnicity, vital_status, last_known_disease_status) %>%
    gather(group_type, group_id, -sample_id)

  sample_groupings <- bind_rows(
    sample_groupings_nci,
    sample_groupings_clin
  ) %>%
    filter(!is.na(group_id), !tolower(group_id) %in% c("not reported", "unknown"))
  # relabel_out <- label_clusters(sample_groupings, dimred_out$cluster, arrange_fun = function(df) filter(PPV > .2) %>% arrange(desc(TPR)))
  relabel_out <- bred::label_clusters(sample_groupings, dimred_out$cluster, arrange_fun = function(df) df %>% arrange(desc(F1)))
  relabel_out$labels %>% print(n = 30)

  # look at cluster annotation

  group_type_palette <-
    unique(sample_groupings$group_type) %>%
    {setNames(RColorBrewer::brewer.pal(length(.), "Dark2"), .)}
  pdf(paste0(data_dir, "cluster_labelling.pdf"), width = 12, height = 6)
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
        facet_wrap(~metric, scales = "free_x", nrow = 1)

      print(g)
    }
  }, finally = {
    dev.off()
  })


  clus <- dimred_out$cluster
  tab <- as.matrix(table(
    paste0(clus, " ", relabel_out$labels$name[clus]),
    sample_info %>% slice(match(names(clus), id)) %>% pull(project_id)
  ))
  tab <- sweep(tab, 1, rowSums(tab), "/")
  tab <- tab[order(apply(tab, 1, max), decreasing = TRUE), ]
  tab <- tab[, order(apply(tab, 2, which.max))]
  pheatmap::pheatmap(tab, angle_col = 315, cluster_rows = FALSE, cluster_cols = FALSE, filename = paste0(data_dir, "cluster_heatmap.pdf"), width = 16, height = 8)

  # save dimred and clustering
  write_rds(c(dimred_out, relabel_out, list(sample_groupings)), paste0(data_dir, "dimred_and_clustering.rds"))
}

list2env(read_rds(paste0(data_dir, "dimred_and_clustering.rds")), .GlobalEnv)

# save visualisation
df <- data.frame(
  dimred_fr,
  sample_info,
  cluster = cluster
) %>%
  left_join(labels, by = "cluster")
labeldf <- df %>% group_by(name) %>% summarise(comp_1 = mean(comp_1), comp_2 = mean(comp_2))

g <-
  ggplot(df, aes(comp_1, comp_2, col = name)) +
  geom_point(size = .5) +
  shadowtext::geom_shadowtext(aes(label = name), labeldf, bg.colour = "white", size = 5) +
  dynplot::theme_graph() +
  coord_equal() +
  scale_colour_manual(values = palette) +
  theme(legend.position = "none")
g
ggsave(paste0(data_dir, "plot_fr.pdf"), g, width = 15, height = 13)

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
ggsave(paste0(data_dir, "legend.pdf"), g, width = 15, height = 6)

rm(expr)
gc()

imp_grouped <-
  importance_sc %>%
  mutate(cluster = as.vector(cluster[sample_id])) %>%
  left_join(labels %>% select(cluster, name), by = "cluster") %>%
  group_by(name, interaction_id, regulator, target) %>%
  summarise(
    importance = mean(importance),
    importance_sc = mean(importance_sc),
    effect = effect[[1]]
  ) %>%
  ungroup() %>%
  left_join(labels %>% select(name, col), by = "name") %>%
  select(source = regulator, name, target, everything()) %>%
  arrange(desc(importance))

imp_grouped_f <- imp_grouped %>% group_by(name) %>% slice(1:50) %>% ungroup() %>% filter(importance > .15)
imp_grouped_f %>% group_by(name) %>% summarise(n = n()) %>% arrange(desc(n)) %>% print(n = Inf)
write_tsv(imp_grouped_f, paste0(data_dir, "grouped_interactions.tsv"))

gc()



# sample_ids <- levels(importance_sc$sample_id)
# interaction_ids <- levels(importance_sc$interaction_id)
# use_scaled_imp <- FALSE
# imp_sc_mat <- Matrix::sparseMatrix(
#   i = importance_sc$sample_id %>% as.integer,
#   j = importance_sc$interaction_id %>% as.integer,
#   x = importance_sc[[if (use_scaled_imp) "importance_sc" else "importance"]],
#   dims = c(length(sample_ids), length(interaction_ids)),
#   dimnames = list(sample_ids, interaction_ids)
# )
# impannot <- pbapply::pblapply(seq_len(500), cl = 1, function(rowi) {
#   intid <- importance$interaction_id[rowi]
#   x <- imp_sc_mat[,intid]
#   test <-
#     sample_groupings %>%
#     group_by(group_type, group_id) %>%
#     summarise(
#       samples = list(match(sample_id, rownames(imp_sc_mat))),
#       pct = length(sample_id) / nrow(imp_sc_mat)
#     ) %>%
#     rowwise() %>%
#     filter(pct > .01, pct < .99) %>%
#     mutate(
#       test = list(t.test(
#         x[samples],
#         x[-samples],
#         alternative = "greater"
#       )),
#       statistic = test$statistic,
#       p_value = test$p.value
#     ) %>%
#     ungroup() %>%
#     arrange(desc(statistic))
#   labels <-
#     test %>%
#     slice(1) %>%
#     mutate(interaction_id = intid) %>%
#     select(-samples, -pct, -test) %>%
#     left_join(importance %>% slice(rowi), by = "interaction_id")
#   lst(test, labels)
# })
# implabels <- impannot %>% map_df("labels")
# write_tsv(implabels, paste0(data_dir, "annotated_interactions.tsv"))

nb_term <- nciont %>% filter(name == "Neuroblastoma") %>% pull(id)
nciont %>% slice(match(nb_ids, id)) %>% print(n = 30)
derp <- sample_info %>% filter(!nci_ont %in% nb_ids, project_id == "TARGET-NBL") %>% pull(nci_ont) %>% unique()
nciont %>% slice(match(derp, id))

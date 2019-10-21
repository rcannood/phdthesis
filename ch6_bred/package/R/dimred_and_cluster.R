
dimred_and_cluster <- function(importance_sc, lmds_ndim = 20, knn = 100, use_scaled_imp = TRUE, knn_dist_quantile = .95) {
  sample_ids <- levels(importance_sc$sample_id)
  interaction_ids <- levels(importance_sc$interaction_id)

  imp_sc_mat <- Matrix::sparseMatrix(
    i = importance_sc$sample_id %>% as.integer,
    j = importance_sc$interaction_id %>% as.integer,
    x = importance_sc[[if (use_scaled_imp) "importance_sc" else "importance"]],
    dims = c(length(sample_ids), length(interaction_ids)),
    dimnames = list(sample_ids, interaction_ids)
  )

  dimred_lmds <- lmds::lmds(imp_sc_mat, ndim = lmds_ndim, distance_method = "spearman")
  rm(imp_sc_mat)
  gc()

  # compute knn
  knn <- RANN::nn2(dimred_lmds, k = knn + 1)
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
    mutate(weight = 1 / dist) %>%
    filter(dist < quantile(dist, knn_dist_quantile))

  gr <- igraph::graph_from_data_frame(
    knndf %>% select(i, j, dist),
    vertices = seq_len(nrow(dimred_lmds)),
    directed = FALSE
  )
  cl <- igraph::cluster_louvain(gr)
  cluster <- cl$membership
  names(cluster) <- rownames(dimred_lmds)

  dimred_fr <- igraph::layout_with_fr(gr)
  rownames(dimred_fr) <- rownames(dimred_lmds)
  colnames(dimred_fr) <- paste0("comp_", seq_len(ncol(dimred_fr)))

  lst(
    dimred_lmds,
    knn,
    dimred_fr,
    cluster
  )
}

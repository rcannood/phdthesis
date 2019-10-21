`%s/%` <- function(x, y) ifelse(y == 0, 1, x / y)

label_clusters <- function(sample_groupings, cluster, arrange_fun = function(df) arrange(df, desc(F1))) {
  if (!"group_type" %in% colnames(sample_groupings)) {
    sample_groupings$group_type <- factor("group")
  }
  TOT <- length(unique(sample_groupings$sample_id))

  test <-
    enframe(cluster, "sample_id", "cluster") %>%
    group_by(cluster) %>%
    mutate(NUM_CLUSTER = n()) %>%
    ungroup() %>%
    left_join(sample_groupings, by = "sample_id") %>%
    group_by(group_type, group_id) %>%
    mutate(NUM_GROUP = n()) %>%
    group_by(cluster, group_type, group_id) %>%
    summarise(
      NUM_CLUSTER = NUM_CLUSTER[[1]],
      NUM_GROUP = NUM_GROUP[[1]],
      TP = n()
    ) %>%
    ungroup() %>%
    mutate(
      FN = NUM_GROUP - TP,
      FP = NUM_CLUSTER - TP,
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
    arrange_fun() %>%
    group_by(cluster) %>%
    slice(1) %>%
    ungroup() %>%
    rename(name = group_id) %>%
    group_by(name) %>%
    mutate(name2 = if (n() == 1) name else paste0(name, " #", row_number())) %>%
    ungroup() %>%
    mutate(name = factor(name2)) %>%
    select(-name2) %>%
    arrange(name) %>%
    mutate(col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(n())) %>%
    arrange(cluster)

  palette <- labels %>% select(name, col) %>% deframe()

  lst(test, labels, palette)
}

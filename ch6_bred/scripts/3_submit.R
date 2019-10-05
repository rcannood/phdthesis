library(tidyverse)

list2env(read_rds("derived_files/data.rds"), .GlobalEnv)

expr <- as.matrix(dataset$expression)

# calculate_target_importance <- function(
#   target_ix,
#   expr,
#   samples,
#   regulators,
#   targets,
#   num_trees = 10000,
#   num_variables_per_split = 100,
#   num_samples_per_tree = 250,
#   min_node_size = 10,
#   interaction_importance_filter = .01,
#   sigmoid_mean = mean(expr[expr != 0]),
#   sigmoid_sd = sd(expr[expr != 0]),
#   num_perms = 10
# ) {
#
#   target <- targets[[target_ix]]
#
#   regs <- setdiff(regulators, target)
#
#   target_expr <- scale(expr[,target])[,1]
#
#   # data <- Matrix::Matrix(expr[,union(target, regs)], sparse = TRUE)
#   # data[,1] <- target_expr
#   # colnames(data)[[1]] <- "PREDICT"
#
#   data <- data.frame(
#     PREDICT = target_expr,
#     expr[,regs],
#     check.names = FALSE,
#     stringsAsFactors = FALSE
#   )
#   # st1 <- system.time({
#   rf2 <- ranger::ranger(
#     data = data,
#     dependent.variable.name = "PREDICT",
#     verbose = TRUE,
#     num.threads = 1,
#     importance = "impurity",
#     mtry = num_variables_per_split,
#     num.trees = num_trees,
#     min.node.size = min_node_size,
#     sample.fraction = num_samples_per_tree / nrow(expr),
#     keep.inbag = TRUE
#   )
#   # })
#
#   # st2 <- system.time({
#   # rf <- randomForest::randomForest(
#   #   expr[,regs],
#   #   target_expr,
#   #   mtry = num_variables_per_split,
#   #   ntree = num_trees,
#   #   sampsize = num_samples_per_tree,
#   #   nodesize = min_node_size,
#   #   importance = FALSE,
#   #   localImp = FALSE,
#   #   keep.inbag = TRUE
#   # )
#   # })
#
#   imp <-
#     rf2$variable.importance %>%
#     enframe("feature_id", "importance") %>%
#     mutate(ix = row_number()) %>%
#     # rf$importance %>%
#     # as.data.frame() %>%
#     # tibble::rownames_to_column("feature_id") %>%
#     # tibble::as_tibble() %>%
#     # mutate(ix = row_number()) %>%
#     # rename(importance = `%IncMSE`) %>%
#     arrange(desc(importance))
#
#   impf <- imp %>% filter(importance >= interaction_importance_filter * 100) # if IncNodePurity
#   impf$effect <- NA
#
#   reg_check <- impf$feature_id
#   limp <- matrix(0, nrow = length(samples), ncol = length(reg_check), dimnames = list(samples, reg_check))
#
#   # inbag <- rf$inbag
#   inbag <- do.call(cbind, rf2$inbag.counts)
#
#   shuffles <- map(seq_len(num_perms), ~ sample.int(nrow(expr)))
#   pred <- predict(rf2, expr, predict.all = TRUE, num.threads = 1)
#   # pred_ind <- pred$individual
#   pred_ind <- pred$predictions
#   pred_ind[inbag > 0] <- NA
#   pred_ind_sw <- sweep(pred_ind, 1, target_expr, "-")^2
#   pred_agg <- rowMeans(pred_ind, na.rm = TRUE)
#
#   expr2 <- expr
#   eff_x <- eff_y <- rep(NA, nrow(expr) * num_perms)
#
#   for (j in seq_along(reg_check)) {
#     regj <- reg_check[[j]]
#     cat("Running permutations for regulator ", j, " / ", length(reg_check), ": ", regj, "\n", sep = "")
#
#     expr_regj <- expr[,regj]
#
#     for (i in seq_len(num_perms)) {
#       ix <- shuffles[[i]]
#
#       expr2[,regj] <- expr_regj[ix]
#
#       pred2 <- predict(rf2, expr2, predict.all = TRUE, num.threads = 1)
#       # pred2_ind <- pred2$individual
#       pred2_ind <- pred2$predictions
#       pred2_ind[inbag > 0] <- NA
#
#       inc_se <-
#         sweep(pred2_ind, 1, target_expr, "-")^2 - pred_ind_sw
#
#       # calculate importance
#       imp <- rowMeans(inc_se, na.rm = TRUE)
#
#       limp[, j] <- limp[, j] + imp / num_perms
#
#       # determine effect
#       effix <- seq_len(nrow(expr)) + nrow(expr) * (i - 1)
#       eff_x[effix] <- expr2[,regj] - expr[,regj]
#       eff_y[effix] <- rowMeans(pred2_ind, na.rm = TRUE) - pred_agg
#     }
#     impf$effect[[j]] <- cor(eff_x, eff_y)
#
#     # reset expr
#     expr2[,regj] <- expr_regj
#   }
#
#   # we've now computed the permutation importance
#   impf$importance <- colMeans(limp)
#
#   # downscale importance if regulator is not expressed
#   # ... it can't be regulating anything if it is not expressed ...
#   expr_reg_sc <- stats::pnorm(expr[,reg_check], mean = sigmoid_mean, sd = sigmoid_sd)
#   limp_sc <- limp * expr_reg_sc
#
#   limp_df <- limp_sc %>%
#     reshape2::melt(varnames = c("cell_id", "regulator"), value.name = "importance") %>%
#     as_tibble()
#
#   lst(
#     importance = impf %>% transmute(
#       regulator = factor(feature_id, levels = regulators),
#       target = factor(target, levels = targets),
#       importance,
#       effect
#     ),
#     importance_sc = limp_df %>% transmute(
#       cell_id = factor(as.character(cell_id), levels = samples),
#       regulator = factor(as.character(regulator), levels = regulators),
#       target = factor(target, levels = targets),
#       importance = importance
#     )
#   )
# }
#
# x <- 1
# qsub_handle <- qsub::qsub_lapply(
#   X = seq_along(targets),
#   qsub_config = qsub::override_qsub_config(
#     # qsub_config = qsub_config,
#     memory = "10G",
#     name = "bred",
#     wait = FALSE,
#     stop_on_error = FALSE,
#     remove_tmp_folder = FALSE,
#     # batch_tasks = round(length(targets) / 1000) # split work into ~ 1000 tasks
#   ),
#   qsub_packages = c("dynutils", "dplyr", "purrr", "magrittr", "tibble"),
#   qsub_environment = c("x"),
#   FUN = calculate_target_importance,
#   # pass data and other parameters
#   expr = expr,
#   samples = samples,
#   regulators = regulators,
#   targets = targets,
#   num_trees = num_trees,
#   num_variables_per_split = num_variables_per_split,
#   num_samples_per_tree = num_samples_per_tree,
#   min_node_size = min_node_size,
#   interaction_importance_filter = interaction_importance_filter,
#   sigmoid_mean = sigmoid_mean,
#   sigmoid_sd = sigmoid_sd,
#   num_perms = num_perms
# )
# write_rds(qsub_handle, "derived_files/qsub_handle.rds")

# qsub_handle <- read_rds("derived_files/qsub_handle.rds")
#
# grn <- qsub::qsub_retrieve(
#   qsub_handle,
#   wait = "just_do_it"
# )
# write_rds(grn, "derived_files/grn.rds")

grn <- read_rds("derived_files/grn.rds")

# remove unfinished executions
grn[map_lgl(grn, ~ length(.) == 1 && is.na(.))] <- NULL

# return combined results
importance <- grn %>% map_df("importance") %>% filter(importance > .1) %>% mutate(i = row_number(), name = paste0(regulator, "->", target))
write_tsv(importance, "~/importance.tsv")
importance_sc <- grn %>%
  map_df("importance_sc") %>%
  filter(importance > .1) %>%
  inner_join(importance %>% select(-importance, -name), by = c("regulator", "target"))

imp_sc_mat <- Matrix::sparseMatrix(
  i = importance_sc$cell_id %>% as.integer,
  j = importance_sc$i,
  x = importance_sc$importance,
  dimnames = list(samples, importance$name)
)

dimred <- dyndimred::dimred_landmark_mds(imp_sc_mat)
g1 <- dynplot::plot_dimred(model, grouping = dataset$grouping) + ggtitle("orig dimred")
g2 <- dynplot::plot_dimred(model, dimred = dimred, grouping = dataset$grouping) + ggtitle("new dimred")
patchwork::wrap_plots(g1, g2, nrow = 1)
dynplot::plot_heatmap(model, imp_sc_mat, features_oi = 100, grouping = dataset$grouping)

dimred2 <- dyndimred::dimred_landmark_mds(Matrix::t(imp_sc_mat), ndim = 10)
cl <- kmeans(dimred2, centers = 8)

rgl::plot3d(dimred2, col = cl$cluster, size = 10)

plot_mat <- imp_sc_mat %>% as.matrix %>% dynutils::scale_quantile() %>% t()
lin <- dynplot:::linearise_cells(model)$progressions %>% arrange(edge_id, cumpercentage)
gaps_col <- which(diff(as.integer(lin$edge_id)) != 0)
pheatmap::pheatmap(
  plot_mat[order(cl$cluster),lin$cell_id],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  gaps_row = which(diff(sort(cl$cluster)) != 0),
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = dataset$grouping %>% enframe("name", "annoation") %>% as.data.frame() %>% column_to_rownames("name")
)

grouping <- model$progressions %>%
  mutate(group = case_when(
    from == "3" & percentage < .33 ~ "MEF",
    from == "3" & percentage < .66 ~ "early-induced",
    from == "3" ~ "induced",
    percentage < .33 ~ "induced",
    to == "4" & percentage < .66 ~ "early-neuron",
    to == "5" & percentage < .66 ~ "early-myocyte",
    to == "4" ~ "neuron",
    to == "5" ~ "myocyte"
  ))

imp_gr <-
  importance_sc %>%
  left_join(grouping %>% select(cell_id, group), by = "cell_id") %>%
  group_by(group, regulator, target) %>%
  summarise_at(c("importance", "effect"), mean) %>%
  ungroup()


imp_gr %>%
  filter(importance > 1) %>%
  write_tsv("~/interactions.tsv")

imp_gr %>%
  group_by(group) %>%
  top_n(50, importance) %>%
  write_tsv("~/interactions_top50.tsv")

plot_dimred(model, grouping = grouping %>% select(cell_id, group) %>% deframe())

# library(tidygraph)
# library(ggraph)
# graph <- as_tbl_graph(importance %>% filter(importance > .1))
#
# # plot using ggraph
# arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(7, "mm"))
# arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(7, "mm"))
# ggraph(graph, layout = 'kk') +
#   geom_edge_fan(aes(colour = effect, width = (importance), filter = effect >= .1), arrow = arrow_up) +
#   geom_edge_fan(aes(colour = effect, width = (importance), filter = effect <= -.1), arrow = arrow_down) +
#   geom_edge_fan(aes(colour = effect, width = (importance), filter = effect > -.1 & effect < .1)) +
#   geom_node_circle(aes(r = .2), fill = "white") +
#   geom_node_text(aes(label = name)) +
#   scale_edge_width_continuous(range = c(.5, 3)) +
#   scale_edge_colour_distiller(palette = "RdBu", direction = 1) +
#   coord_equal() +
#   theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

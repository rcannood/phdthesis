library(dlstats)
library(tidyverse)

cran <- c(
  "babelwhale", "diffusionMap", "dyndimred", "dynparam", "dynutils", "dynwrap",
  "GillespieSSA", "GillespieSSA2", "incgraph", "lmds", "princurve", "proxyC",
  "qsub", "SCORPIUS", "badger", "devtools", "ggrepel", "mlr", "ParamHelpers", "ranger", "remotes", "rlang", "shadowtext", "tidyverse"
)
bioc <- c("ClusterSignificance", "monocle", "slingshot", "splatter")

dl_cran <- dlstats::cran_stats(cran)
dl_bioc <- dlstats::bioc_stats(bioc)

dls <- bind_rows(dl_cran, dl_bioc %>% select(start, end, package, downloads = Nb_of_downloads)) %>%
  mutate(role = ifelse(package %in% c("devtools", "mlr", "ParamHelpers", "remotes", "rlang", "shadowtext", "ggrepel", "tidyverse", bioc), "Contributor", "Author"))

min_date <- "2019-10-01"
num_days <- as.integer(Sys.Date() - as.Date(min_date))
summ <- dls %>% filter(min_date <= start) %>% group_by(package) %>% summarise(downloads = sum(downloads) / num_days * 365) %>% slice(match(cran, package))
summ %>% print(n = 50)

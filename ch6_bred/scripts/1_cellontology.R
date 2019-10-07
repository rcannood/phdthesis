library(tidyverse)
library(ontologyIndex)


cellontf <- "derived_files/cl.obo"

if (!file.exists(cellontf)) download.file("https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl.obo", cellontf)

cell_ont <-
  get_ontology(cellontf) %>%
  with(tibble(id, name, parents, children, ancestors, obsolete)) %>%
  filter(!obsolete, grepl("CL:", id))

gr <- cell_ont %>% select(from = id, to = children) %>% unnest(to) %>%
  igraph::graph_from_data_frame()

dis <- igraph::distances(gr, v = "CL:0000000")
depth <- dis[1,] %>% enframe("id", "depth")

cell_ont <- cell_ont %>% left_join(depth, by = "id")

cell_ont_up <-
  cell_ont %>%
  select(cell_ontology_id = id, upstream = ancestors) %>%
  unnest(upstream) %>%
  filter(cell_ontology_id != upstream)

# cids <- unique(sample_info$cell_ontology_id)
cids <- c(
  "CL:0000763", "CL:0000583", "CL:0000236", "CL:0000623", "CL:0000084", "CL:0000875", "CL:0000738", "CL:0000860", "CL:0000097", "CL:0000576", "CL:0002191", "CL:0000559", "CL:0000094", "CL:0000765",
  "CL:0008001", "CL:0000547", "CL:0002048", "CL:0000767", "CL:0000235", "CL:0002046", "CL:0000816", "CL:0002045", "CL:0000453", "CL:0000451", "CL:0000081", "CL:0002420", "CL:0000894"
)
ciddis <- igraph::distances(gr, v = cids, to = cids)
ciddr <- stats::cmdscale(as.dist(ciddis), k = 10)

out <- "derived_files/cell_ontology.rds"
write_rds(lst(cell_ont, cell_ont_up, ciddr), out, compress = "gz")

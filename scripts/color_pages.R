library(tidyverse)

col <- system("gs -o - -sDEVICE=inkcov 'boek_v2.1_20191211.pdf'", intern = TRUE)

pagix <- grep("^Page [0-9]*$", col)
colix <- grep("CMYK OK$", col)

pages <- col[pagix] %>% gsub("Page ", "", .) %>% as.integer()
colors <- col[colix] %>% gsub("^ *(.*) *CMYK OK$", "\\1", .) %>% strsplit("  *") %>% map(as.numeric) %>% do.call(rbind, .)
colnames(colors) <- c("C", "M", "Y", "K")

coldf <- data.frame(
  page = pages,
  colors
)


colpag <- coldf %>% filter(C > 0 | M > 0 | Y > 0)
pages <- colpag$page
paste(pages, collapse = ",")


length(pages)
nrow(coldf) - length(pages)

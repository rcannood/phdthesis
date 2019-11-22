library(tidyverse)
library(scales)

set.seed(1)
df <- tibble(
  x = round(10^runif(6, .5, )),
  n = paste0("Molecule ", seq_along(x)) %>% forcats::fct_rev()
)

ggplot(df) +
  geom_bar(aes(n, x), stat = "identity") +
  theme_classic() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = 10^c(0, 5)) +
  labs(x = NULL, y = "Abundance") +
  coord_flip()


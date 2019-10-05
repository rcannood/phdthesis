library(tidyverse)

library(hdf5r)
fil <- H5File$new("derived_files/scanpy_processed_matrix_both_no_cc.h5ad")

hdf5r:::print_listing(fil, 20)


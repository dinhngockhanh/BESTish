suppressPackageStartupMessages({
    library(BESTish)
    library(dplyr)
    library(tidyr)
    library(readxl)
})
path <- system.file("extdata", "41586_2022_4785_MOESM5_ESM.xlsx", package = "BESTish")
fabre <- readxl::read_xlsx(path, sheet = 1)

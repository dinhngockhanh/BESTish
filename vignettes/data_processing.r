suppressPackageStartupMessages({
    devtools::load_all() ############################################### REPLACE WITH library(BESTish)
    library(dplyr)
    library(tidyr)
    library(readxl)
})
#----------------------Load raw data from Watson et al. and Fabre et al.
path <- here::here("data", "41586_2022_4785_MOESM5_ESM.xlsx") ########## REPLACE WITH path <- system.file("data", "41586_2022_4785_MOESM5_ESM.xlsx", package = "BESTish")
fabre <- readxl::read_xlsx(path, sheet = 1)
fabre <- as.data.frame(fabre, stringsAsFactors = FALSE)
colnames(fabre) <- c(
    "Sample_ID", "Sex", "Study_phase", "Age", "Chromosome", "Start", "End", "WT", "MT", "VAF", "Gene_ID", "Protein_change_ID", "Effect"
)
path <- here::here("data", "all_studies_trimmed_all_genes.csv") ######## REPLACE WITH path <- system.file("data", "all_studies_trimmed_all_genes.csv", package = "BESTish")
watson <- read.csv(path, stringsAsFactors = FALSE)
colnames(watson) <- c(
    "VAF", "Age", "Protein_change_ID", "Gene_ID", "Study"
)
watson$Protein_change_ID <- paste0("p.", watson$Protein_change_ID)

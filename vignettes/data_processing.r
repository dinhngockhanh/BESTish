suppressPackageStartupMessages({
    library(BESTish)
    library(dplyr)
    library(tidyr)
    library(readxl)
})
#----------------------Load raw data from Fabre et al. and Watson et al.
path <- system.file("extdata", "41586_2022_4785_MOESM5_ESM.xlsx", package = "BESTish")
fabre <- read_xlsx(path, sheet = 1)
fabre <- as.data.frame(fabre, stringsAsFactors = FALSE)
colnames(fabre) <- c(
    "Sample_ID", "Sex", "Study_phase", "Age", "Chromosome", "Start", "End", "WT", "MT", "VAF", "Gene_ID", "Protein_change_ID", "Effect"
)
path <- system.file("extdata", "all_studies_trimmed_all_genes.csv", package = "BESTish")
watson <- read.csv(path, stringsAsFactors = FALSE)
colnames(watson) <- c(
    "VAF", "Age", "Protein_change_ID", "Gene_ID", "Study"
)
watson$Protein_change_ID <- paste0("p.", watson$Protein_change_ID)
#-----------------Determine potential CH drivers to analyze with BESTish
#----------------as mutations present at >=3 time points in >=3 patients
#---------------------------------in longitudinal data from Fabre et al.
potential_CH_drivers <- fabre %>%
    dplyr::group_by(Sample_ID, Gene_ID, Protein_change_ID) %>%
    dplyr::summarise(n_ages = dplyr::n_distinct(Age), .groups = "drop") %>%
    dplyr::filter(n_ages >= 3) %>%
    dplyr::group_by(Gene_ID, Protein_change_ID) %>%
    dplyr::summarise(n_sard = dplyr::n_distinct(Sample_ID), .groups = "drop") %>%
    dplyr::filter(n_sard >= 3) %>%
    dplyr::select(Gene_ID, Protein_change_ID) %>%
    dplyr::distinct()
fabre_inference_list <- fabre %>%
    dplyr::semi_join(potential_CH_drivers, by = c("Gene_ID", "Protein_change_ID")) %>%
    dplyr::distinct(Protein_change_ID, Sample_ID, Gene_ID)
write.csv(fabre, file = "fabre.csv", row.names = FALSE)
write.csv(fabre_inference_list, file = "fabre_inference_list.csv", row.names = FALSE)
#---------Determine which potential CH drivers have adequate cohort data
#---------------------(>= 8 patients) in Coombs et al. or Watson et al.
watson_inference_list <- watson %>%
    dplyr::semi_join(
        fabre_inference_list %>% dplyr::select(Protein_change_ID, Gene_ID) %>% dplyr::distinct(),
        by = c("Protein_change_ID", "Gene_ID")
    ) %>%
    dplyr::filter(Study %in% c("Coombs2017", "McKerrel2015")) %>%
    dplyr::group_by(Protein_change_ID, Gene_ID, Study) %>%
    dplyr::filter(dplyr::n() >= 8) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Protein_change_ID, Study, Gene_ID)
write.csv(watson, file = "watson.csv", row.names = FALSE)
write.csv(watson_inference_list, file = "watson_inference_list.csv", row.names = FALSE)

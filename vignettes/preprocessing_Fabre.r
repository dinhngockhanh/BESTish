suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readxl)
})
setwd("/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")
path <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/data/41586_2022_4785_MOESM5_ESM.xlsx"
fabre <- readxl::read_xlsx(path, sheet = 1)
fabre <- as.data.frame(fabre, stringsAsFactors = FALSE)
colnames(fabre) <- c(
    "Sample_ID", "Sex", "Study_phase", "Age", "Chromosome", "Start", "End", "WT", "MT", "VAF", "Gene_ID", "Protein_change_ID", "Effect"
)
path <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/data/all_studies_trimmed_all_genes.csv"
watson <- read.csv(path, stringsAsFactors = FALSE)
colnames(watson) <- c(
    "VAF", "Age", "Protein_change_ID", "Gene_ID", "Study"
)
watson$Protein_change_ID <- paste0("p.", watson$Protein_change_ID)
#---Keep only unique variants present at >=3 age points in >=3 patients
variants_kept <- fabre %>%
    dplyr::group_by(Sample_ID, Gene_ID, Protein_change_ID) %>%
    dplyr::summarise(n_ages = dplyr::n_distinct(Age), .groups = "drop") %>%
    dplyr::filter(n_ages >= 3) %>%
    dplyr::group_by(Gene_ID, Protein_change_ID) %>%
    dplyr::summarise(n_sard = dplyr::n_distinct(Sample_ID), .groups = "drop") %>%
    dplyr::filter(n_sard >= 3) %>%
    dplyr::select(Gene_ID, Protein_change_ID) %>%
    dplyr::distinct()
fabre_filtered <- fabre %>%
    dplyr::semi_join(variants_kept, by = c("Gene_ID", "Protein_change_ID"))
fabre_inference_list <- fabre %>%
    dplyr::distinct(Protein_change_ID, Sample_ID, Gene_ID)
fabre_inference_list_filtered <- fabre_filtered %>%
    dplyr::distinct(Protein_change_ID, Sample_ID, Gene_ID)
write.csv(fabre, file = "fabre_master.csv", row.names = FALSE)
write.csv(fabre_filtered, file = "fabre_master_filtered.csv", row.names = FALSE)
write.csv(fabre_inference_list, file = "fabre_inference_list.csv", row.names = FALSE)
write.csv(fabre_inference_list_filtered, file = "fabre_inference_list_filtered.csv", row.names = FALSE)
cat("Unique Protein_change_ID x Sample_ID pairs in Fabre:               ", nrow(fabre_inference_list), "\n")
cat("Unique Protein_change_ID x Sample_ID pairs in filtered Fabre:      ", nrow(fabre_inference_list_filtered), "\n")
#---Keep only cohort data for variants in filtered list and >= 10 patients per study
watson_filtered <- watson %>%
    dplyr::semi_join(
        fabre_inference_list_filtered %>% dplyr::select(Protein_change_ID, Gene_ID) %>% dplyr::distinct(),
        by = c("Protein_change_ID", "Gene_ID")
    ) %>%
    dplyr::filter(Age != "noagedata") %>%
    dplyr::group_by(Protein_change_ID, Gene_ID, Study) %>%
    dplyr::filter(dplyr::n() >= 8) %>%
    dplyr::ungroup()
watson_inference_list <- watson %>%
    dplyr::distinct(Protein_change_ID, Study, Gene_ID)
watson_inference_list_filtered <- watson_filtered %>%
    dplyr::distinct(Protein_change_ID, Study, Gene_ID)
write.csv(watson, file = "watson_master.csv", row.names = FALSE)
write.csv(watson_filtered, file = "watson_master_filtered.csv", row.names = FALSE)
write.csv(watson_inference_list, file = "watson_inference_list.csv", row.names = FALSE)
write.csv(watson_inference_list_filtered, file = "watson_inference_list_filtered.csv", row.names = FALSE)
cat("Unique Protein_change_ID x Study pairs in Watson:                  ", nrow(watson_inference_list), "\n")
cat("Unique Protein_change_ID x Study pairs in filtered Watson:         ", nrow(watson_inference_list_filtered), "\n")

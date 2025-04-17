library(dplyr)
library(readr)
# BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
# ======================================================================
# #---Input COSMIC mutation data
# cosmic_file <- paste0(wd, "data/COSMIC/Cosmic_NonCodingVariants_v101_GRCh38.vcf")
# cosmic_data <- readVcf(vcf_file, "hg38")
# cosmic_df<- data.frame(
#     LEGACY_ID=info(cosmic_data)$LEGACY_ID,
#     GENE=info(cosmic_data)$GENE
# )
# if(any(is.na(cosmic_df$GENE))) cosmic_df<- cosmic_df[!is.na(cosmic_df$GENE),]
# cosmic_df <- cosmic_df %>%
#     group_by(LEGACY_ID) %>%
#     summarise(GENE = paste(unique(GENE), collapse = ","), .groups = "drop")
#---Input clinical information for all ICGC samples
icgc_data <- read_csv(file = "data/PCAWG/ICGC_sample_information.csv", guess_max = 100000)
files_in_dir <- list.files("data/PCAWG", full.names = TRUE)
matching_files <- files_in_dir[grepl("_all\\.csv$", files_in_dir)]
file_ids <- sub("_all\\.csv$", "", basename(matching_files))
#---Prepare data for every ICGC cancer type
dir.create("vignettes/PCAWG", recursive = TRUE, showWarnings = FALSE)
for (cancer_type in unique(icgc_data$histology_abbreviation)) {
    dir.create(paste0("vignettes/PCAWG_",cancer_type), recursive = TRUE, showWarnings = FALSE)
    cancer_type_data <- icgc_data %>%
        filter(
            histology_abbreviation == cancer_type,
            aliquot_id %in% file_ids
        )
    output_file <- file.path("vignettes", "PCAWG", paste0(cancer_type, "_clinical_information.csv"))
    write_csv(cancer_type_data, output_file)
    matched_ids <- unique(cancer_type_data$aliquot_id)
    for (aliquot in matched_ids) {
        csv_file <- file.path(all_csv_dir, paste0(aliquot, "_all.csv"))
        if (file.exists(csv_file)) {
            all_data <- read_csv(file = csv_file, guess_max = 100000)
            filtered_data <- all_data %>% filter(!is.na(cosmic))
            filtered_data$sample <- aliquot
            filtered_data$gene <- NA
            filtered_data$gene <- cosmic_df$GENE[match(filtered_data$cosmic, cosmic_df$LEGACY_ID)]
            output_file <- file.path("vignettes", paste0("PCAWG_", cancer_type), paste0(aliquot, ".csv"))
            write_csv(filtered_data, output_file)
        }
    }
}

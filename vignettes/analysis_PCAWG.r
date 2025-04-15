library(dplyr)
library(readr)
library(VariantAnnotation)
# ------------------------------------------ Working directory for Keito
wd <- "/Users/keitotaketomi/Documents/DriverSelectionSweep/"
# ------------------------------------------ Working directory for Khanh
wd <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/"
# ======================================================================
cancer_type <- "Ovary-AdenoCA"
#---Input COSMIC mutation data
vcf_file <- paste0("data/COSMIC/Cosmic_NonCodingVariants_Vcf_v101_GRCh38/Cosmic_NonCodingVariants_v101_GRCh38.vcf")
vcf_data <- readVcf(vcf_file, "hg38")


#---Input mutational information from ICGC
icgc_file <- paste0(wd, "/data/PCAWG/ICGC_sample_information.csv")
icgc_data <- read_csv(file = icgc_file, guess_max = 100000)
all_csv_dir <- paste0(wd, "/data/PCAWG")
files_in_dir <- list.files(all_csv_dir, full.names = TRUE)
matching_files <- files_in_dir[grepl("_all\\.csv$", files_in_dir)]
file_ids <- sub("_all\\.csv$", "", basename(matching_files))
cancer_type_data <- icgc_data %>%
    filter(
        histology_abbreviation == cancer_type,
        aliquot_id %in% file_ids
    )
matched_ids <- unique(cancer_type_data$aliquot_id)
for (aliquot in matched_ids) {
    csv_file <- file.path(all_csv_dir, paste0(aliquot, "_all.csv"))
    if (file.exists(csv_file)) {
        all_data <- read_csv(file = csv_file, guess_max = 100000)
        filtered_data <- all_data %>% filter(!is.na(cosmic))
    }
}

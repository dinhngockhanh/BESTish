library(VariantAnnotation)
vcf_file <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/COSMIC/Cosmic_NonCodingVariants_Vcf_v101_GRCh38/Cosmic_NonCodingVariants_v101_GRCh38.vcf"
vcf_data <- readVcf(vcf_file, "hg38")
summary(vcf_data)

legacy_id <- "COSN28762392"

matching_indices <- which(info(vcf_data)$LEGACY_ID == legacy_id)
matching_info <- info(vcf_data)[matching_indices, , drop = FALSE]
matching_row_data <- rowRanges(vcf_data)[matching_indices, ]

library(readxl)
library(dplyr)
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG"
#-------------------------------------------------------Input PCAWG data
specimen_histology_data <- read_excel(paste0(R_data, "/pcawg_specimen_histology_August2016_v9.xlsx"))
sample_data <- read_excel(paste0(R_data, "/pcawg_sample_sheet.xlsx"))


specimen_histology_data <- specimen_histology_data[
    which(specimen_histology_data$histology_abbreviation == "ColoRect-AdenoCA"),
]

combined_data <- specimen_histology_data %>%
    left_join(sample_data, by = "icgc_sample_id") %>%
    filter(library_strategy == "WGS")


for (aliquot_id in combined_data$aliquot_id) {
    print(aliquot_id)
    CNV_path <- paste0(R_data, "/consensus_cnv/consensus.20170119.somatic.cna.annotated/", aliquot_id, ".consensus.20170119.somatic.cna.annotated.txt")
    SNV_path <- paste0(R_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv/", aliquot_id, ".consensus.20160830.somatic.snv_mnv.vcf")
    print(file.exists(CNV_path))
    print(file.exists(SNV_path))
}

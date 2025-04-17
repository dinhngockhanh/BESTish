library(readr)
# ======================================================================
cancer_type <- "Ovary-AdenoCA"
gene <- "TP53"
#---Input clinical information
input_file <- file.path("vignettes", "PCAWG", paste0(cancer_type, "_clinical_information.csv"))
cancer_type_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
inactivation_variant_classifications <- c("Nonsense_Mutation", "Splice_Site", "Missense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del")
cancer_type_data[[paste0("Classification_", gene)]] <- NA
for (row in 1:nrow(cancer_type_data)) {
    aliquot_id <- cancer_type_data$aliquot_id[row]
    input_file <- file.path("vignettes", paste0("PCAWG_", cancer_type), paste0(aliquot_id, ".csv"))
    mutation_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
    gene_mutation_rows <- mutation_data[grepl(paste0("\\b", gene, "\\b"), mutation_data$gene, ignore.case = TRUE), ]
    gene_mutation_rows <- gene_mutation_rows[which(gene_mutation_rows$Variant_Classification %in% inactivation_variant_classifications), ]
    if (nrow(gene_mutation_rows) == 0) {
        sample_classification <- 0
    } else {
        gene_mutation_rows$classification <- ifelse(gene_mutation_rows$minor_cn == 0, 2, 1)
        sample_classification <- max(gene_mutation_rows$classification)
    }
    cancer_type_data[[paste0("Classification_", gene)]][row] <- sample_classification
}
write_csv(cancer_type_data, file.path("vignettes", paste0(cancer_type, "_", gene, ".csv")))

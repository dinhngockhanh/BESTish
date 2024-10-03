library(dplyr)
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/MSK-IMPACT [2017]"
#--------------------------------------------------Input MSK-IMPACT data
clinical_sample_data <- read.csv(paste0(R_data, "/data_clinical_sample_standardized.csv"), header = TRUE)
clinical_patient_data <- read.csv(paste0(R_data, "/data_clinical_patient_standardized.csv"), header = TRUE)
mutation_data <- read.csv(paste0(R_data, "/data_mutations_standardized.csv"), header = TRUE)
copy_number_data <- read.csv(paste0(R_data, "/data_cna_standardized.csv"), header = TRUE)
colnames(copy_number_data) <- gsub("\\.", "-", colnames(copy_number_data))
#-------------------------------Combine clinical sample and patient data
clinical_data <- clinical_sample_data %>%
    left_join(clinical_patient_data, by = "PATIENT_ID")
#-----------------------------------Find sample purity for each mutation
mutation_data <- mutation_data %>%
    left_join(clinical_data %>% select(SAMPLE_ID, TUMOR_PURITY),
        by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")
    ) %>%
    mutate(Purity = TUMOR_PURITY) %>%
    select(-TUMOR_PURITY)
#---------------------------------------Add copy number to mutation data
mutation_data$Copy_Number <- NA
for (row in 1:nrow(mutation_data)) {
    Hugo_Symbol <- mutation_data$Hugo_Symbol[row]
    if (!Hugo_Symbol %in% copy_number_data$Hugo_Symbol) next
    Tumor_Sample_Barcode <- mutation_data$Tumor_Sample_Barcode[row]
    Copy_Number <- copy_number_data[[Tumor_Sample_Barcode]][which(copy_number_data$Hugo_Symbol == Hugo_Symbol)]
    mutation_data$Copy_Number[row] <- Copy_Number
}
#-------------------------------------------------Compute mutational VAF
mutation_data$VAF <- mutation_data$t_alt_count / (mutation_data$t_ref_count + mutation_data$t_alt_count)
#-------------------------------------------------Compute mutational CCF
mutation_data$CCF <- 2 * 100 * mutation_data$VAF / mutation_data$Purity





########################################################################
######################################### ANALYSIS FOR COLORECTAL CANCER
########################################################################
# cancer_type <- "Colorectal Cancer"
# clinical_data <- clinical_data[which(clinical_data$CANCER_TYPE == cancer_type), ]
# sample_IDs <- clinical_data$SAMPLE_ID
# for (sample_ID in sample_IDs) {
#     sample_mutation <- mutation_data[mutation_data$Tumor_Sample_Barcode == sample_ID, ]
#     if ("TP53" %in% sample_mutation$Hugo_Symbol & "APC" %in% sample_mutation$Hugo_Symbol & "KRAS" %in% sample_mutation$Hugo_Symbol) {
#         print(sample_mutation[sample_mutation$Hugo_Symbol %in% c("TP53", "APC", "KRAS"), c("Hugo_Symbol", "VAF", "Variant_Classification", "Copy_Number")])
#     }
# }

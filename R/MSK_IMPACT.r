library(dplyr)
library(ggplot2)
R_data <- "data"
R_workplace <- "R"
#--------------------------------------------------Input MSK-IMPACT data
clinical_sample_data <- read.csv(paste0(R_data, "/data_clinical_sample_standardized.csv"), header = TRUE)
clinical_patient_data <- read.csv(paste0(R_data, "/data_clinical_patient_standardized.csv"), header = TRUE)
mutation_data <- read.csv(paste0(R_data, "/data_mutations_standardized.csv"), header = TRUE)
gene_copy_number_data <- read.csv(paste0(R_data, "/data_cna_standardized.csv"), header = TRUE)
copy_number_data <- read.delim(paste0(R_data, "/msk_impact_2017_segments.seg"), header = TRUE)
colnames(gene_copy_number_data) <- gsub("\\.", "-", colnames(gene_copy_number_data))
#-------------------------------Combine clinical sample and patient data
clinical_data <- clinical_sample_data %>%
    left_join(clinical_patient_data, by = "PATIENT_ID")
copy_number_data <- copy_number_data %>%
    left_join(clinical_data %>% select(SAMPLE_ID, SEX), by = c("ID" = "SAMPLE_ID"))
#----------------------------Limit to study of primary Colorectal Cancer
cancer_type <- "Colorectal Cancer"
clinical_data <- clinical_data %>%
    filter(SAMPLE_TYPE %in% c("Primary"), CANCER_TYPE == cancer_type)
mutation_data <- mutation_data[mutation_data$Tumor_Sample_Barcode %in% clinical_data$SAMPLE_ID, ]
copy_number_data <- copy_number_data[copy_number_data$ID %in% clinical_data$SAMPLE_ID, ]
#----------------------------------Get segmented CN data for each sample
copy_number_data <- copy_number_data %>%
    mutate(width = loc.end - loc.start + 1) %>%
    mutate(normal_cn = ifelse(SEX == "Male" & chrom %in% c("X", "Y"), 1, 2)) %>%
    mutate(tumor_cn = round(2^seg.mean * normal_cn))
#-------------------------Get Fraction of Genome Altered for each sample
clinical_data$FGA <- 0
for (i in 1:nrow(clinical_data)) {
    sample_id <- clinical_data$SAMPLE_ID[i]
    filtered_data <- copy_number_data %>%
        filter(ID == sample_id)
    altered_width_sum <- filtered_data %>%
        filter(normal_cn != tumor_cn) %>%
        summarize(total_width = sum(width, na.rm = TRUE)) %>%
        pull(total_width)
    total_width_sum <- filtered_data %>%
        summarize(total_width = sum(width, na.rm = TRUE)) %>%
        pull(total_width)
    clinical_data$FGA[i] <- altered_width_sum / total_width_sum
}
#   Plot distribution of mutation burden vs CNA burden in each sample
png(paste0(R_workplace, "/TMB_vs_FGA.png"), width = 800, height = 600)
ggplot(clinical_data, aes(x = TMB_NONSYNONYMOUS, y = FGA)) +
    geom_point() +
    labs(
        title = NULL,
        x = "Nonsynonymous Tumor Mutation Burden",
        y = "Fraction of Genome Altered"
    ) +
    theme_minimal()
dev.off()
#-----------------------------------Find sample purity for each mutation
mutation_data <- mutation_data %>%
    left_join(clinical_data %>% select(SAMPLE_ID, TUMOR_PURITY),
        by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")
    ) %>%
    mutate(Purity = TUMOR_PURITY) %>%
    select(-TUMOR_PURITY)
#---------------------------------------Add copy number to mutation data
mutation_data$normal_cn <- NA
mutation_data$tumor_cn <- NA
for (row in 1:nrow(mutation_data)) {
    matching_row <- copy_number_data %>%
        filter(
            ID == mutation_data$Tumor_Sample_Barcode[row],
            chrom == mutation_data$Chromosome[row],
            loc.start <= mutation_data$Start_Position[row],
            loc.end > mutation_data$Start_Position[row]
        )
    if (nrow(matching_row) > 0) {
        mutation_data$normal_cn[row] <- matching_row$normal_cn
        mutation_data$tumor_cn[row] <- matching_row$tumor_cn
    }
}




#-------------------------------------
mutation_data <- mutation_data %>%
    filter(Hugo_Symbol %in% c("KRAS", "APC", "TP53"))
mutation_data <- mutation_data %>%
    filter(!is.na(Purity) & Purity >= 50 & !is.na(normal_cn) & !is.na(tumor_cn) & tumor_cn <= 2)
#-------------------------------------------------Compute mutational VAF
mutation_data$t_alt_count <- pmax(mutation_data$t_alt_count, 0)
mutation_data$t_ref_count <- pmax(mutation_data$t_ref_count, 0)
mutation_data$t_tot_count <- mutation_data$t_alt_count + mutation_data$t_ref_count
mutation_data$VAF <- mutation_data$t_alt_count / mutation_data$t_tot_count
#-------------------------------------------------Compute mutational CCF
mutation_data$purity <- mutation_data$Purity / 100
mutation_data$multiplicity <- pmax(1, round(mutation_data$VAF / (mutation_data$purity) * (mutation_data$purity * mutation_data$tumor_cn + (1 - mutation_data$purity) * mutation_data$normal_cn)))
mutation_data$CCF <- mutation_data$VAF / (mutation_data$purity * mutation_data$multiplicity) * (mutation_data$purity * mutation_data$tumor_cn + (1 - mutation_data$purity) * mutation_data$normal_cn)
cat("\nBEFORE FILTERING:--------------------------------------------------\n")
cat("Max VAF:           ", max(mutation_data$VAF), "\n")
cat("Min VAF:           ", min(mutation_data$VAF), "\n")
cat("Max CCF:           ", max(mutation_data$CCF), "\n")
cat("Min CCF:           ", min(mutation_data$CCF), "\n")
cat("% with CCF>1.5:    ", length(which(mutation_data$CCF > 1.5)) / nrow(mutation_data), "\n")
# mutation_data <- mutation_data %>%
#     filter(CCF <= 1.5) %>%
#     mutate(CCF = ifelse(CCF > 1, 1, CCF))




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

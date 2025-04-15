library(dplyr)   
library(readr)   

icgc_file <- "/Users/keitotaketomi/Documents/DriverSelectionSweep/data/PCAWG/ICGC_sample_information.csv"

all_csv_dir <- "/Users/keitotaketomi/Documents/DriverSelectionSweep/data/PCAWG"

icgc_data <- read_csv(file = icgc_file, guess_max = 100000)

print(colnames(icgc_data))

files_in_dir <- list.files(all_csv_dir, full.names = TRUE)

matching_files <- files_in_dir[grepl("_all\\.csv$", files_in_dir)]
message("Files matching '_all.csv':")
print(matching_files)

file_ids <- sub("_all\\.csv$", "", basename(matching_files))
message("Extracted Sample IDs from filenames:")
print(file_ids)

ovary_data <- icgc_data %>% 
  filter(histology_abbreviation == "Ovary-AdenoCA",
         aliquot_id %in% file_ids)

if (nrow(ovary_data) == 0) {
  stop("No Ovary-AdenoCA records with matching _all.csv files found. Check sample IDs and file naming.")
}

matched_ids <- unique(ovary_data$aliquot_id)
message("Matched aliquot IDs to process:")
print(matched_ids)

for (aliquot in matched_ids) {
  # Construct the full path for this aliquot's _all.csv file
  csv_file <- file.path(all_csv_dir, paste0(aliquot, "_all.csv"))
  
  if (file.exists(csv_file)) {
    message("Processing file: ", csv_file)
    
    # Read the _all.csv file 
    all_data <- read_csv(file = csv_file, guess_max = 100000)
    
    # Filter: keep rows where 'cosmic' is not NA
    filtered_data <- all_data %>% filter(!is.na(cosmic))
    
    # Write the filtered data to a new CSV file in the same directory
    output_file <- file.path(all_csv_dir, paste0(aliquot, "_all_filtered.csv"))
    write_csv(filtered_data, output_file)
    
    message("Filtered data written to: ", output_file)
  } else {
    warning("File not found for aliquot: ", aliquot, " (Expected: ", csv_file, ")")
  }
}

message("Processing complete.")

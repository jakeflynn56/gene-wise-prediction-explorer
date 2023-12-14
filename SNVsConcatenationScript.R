# Load necessary library
library(dplyr)

# List of .gz file names (replace these with your actual file names)
file_names <- c("BayesDel_170824_noAF_chr1.gz", "BayesDel_170824_noAF_chr2.gz", "BayesDel_170824_noAF_chr3.gz", "BayesDel_170824_noAF_chr4.gz", "BayesDel_170824_noAF_chr5.gz", "BayesDel_170824_noAF_chr6.gz", "BayesDel_170824_noAF_chr7.gz", "BayesDel_170824_noAF_chr8.gz", "BayesDel_170824_noAF_chr9.gz", "BayesDel_170824_noAF_chr10.gz", "BayesDel_170824_noAF_chr11.gz", "BayesDel_170824_noAF_chr12.gz", "BayesDel_170824_noAF_chr13.gz", "BayesDel_170824_noAF_chr14.gz", "BayesDel_170824_noAF_chr15.gz", "BayesDel_170824_noAF_chr16.gz", "BayesDel_170824_noAF_chr17.gz", "BayesDel_170824_noAF_chr18.gz", "BayesDel_170824_noAF_chr19.gz", "BayesDel_170824_noAF_chr20.gz", "BayesDel_170824_noAF_chr21.gz", "BayesDel_170824_noAF_chr22.gz", "BayesDel_170824_noAF_chrM.gz", "BayesDel_170824_noAF_chrX.gz", "BayesDel_170824_noAF_chrY.gz")

# Initialize an empty data frame to store all data
all_data <- data.frame()

# Loop through each file, read the data, and concatenate
for(file_name in file_names) {
  # Read the gz file
  data <- read.csv(gzfile(file_name), header = TRUE, sep = "\t") # Set header=FALSE if there is no header
  data$X.Chr <- as.character(data$X.Chr)  # Convert X.Chr to character
  print(paste("Processing file", file_name))

  # Concatenate the data
  all_data <- bind_rows(all_data, data)
}

# Write the combined data to a new CSV file
write.csv(all_data, "combined_chromosomes_bd.csv", row.names = FALSE)

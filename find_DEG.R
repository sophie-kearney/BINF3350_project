######
# LOAD PACKAGES
######
library(dplyr)
library(janitor)

######
# READ IN DATA
######
luad_data <- read.table(file="LUAD_Differential_Gene_Expression_Table.txt", 
                   header=TRUE, sep="\t")
luad_data <- luad_data %>%
  clean_names()
luad_p_val <- luad_data[luad_data$fdr_adjusted_p_value<=.05,]

brca_data <- read.table(file="BRCA_Differential_Gene_Expression_Table.txt", 
                   header=TRUE, sep="\t")
brca_data <- brca_data %>%
  clean_names()
brca_p_val <- brca_data[brca_data$fdr_adjusted_p_value<=.05,]

kirc_data <- read.table(file="KIRC_Differential_Gene_Expression_Table.txt", 
                   header=TRUE, sep="\t")
kirc_data <- kirc_data %>%
  clean_names()
kirc_p_val <- kirc_data[kirc_data$fdr_adjusted_p_value<=.05,]

paad_data <- read.table(file="PAAD_Differential_Gene_Expression_Table.txt", 
                   header=TRUE, sep="\t")
paad_data <- paad_data %>%
  clean_names()
paad_p_val <- paad_data[paad_data$fdr_adjusted_p_value<=.05,]


######
# FIND COMMON GENES 
######

# remove cancer type column
paad <- paad_p_val[-1]
kirc <- kirc_p_val[-1]
brca <- brca_p_val[-1]
luad <- luad_p_val[-1]

# create unique column names for each cancer
colnames(paad) <- c("ncbi_gene_id","PAAD_fdr_adjusted_p_value", "PAAD_cancer_sample_med",
                    "PAAD_normal_sample_med", "PAAD_log2_fold_change", "PAAD_p_value", 
                    "gene_symbol")
colnames(kirc) <- c("ncbi_gene_id","KIRC_fdr_adjusted_p_value", "KIRC_cancer_sample_med",
                    "KIRC_normal_sample_med", "KIRC_log2_fold_change", "KIRC_p_value", 
                    "gene_symbol")
colnames(brca) <- c("ncbi_gene_id","BRCA_fdr_adjusted_p_value", "BRCA_cancer_sample_med",
                    "BRCA_normal_sample_med", "BRCA_log2_fold_change", "BRCA_p_value", 
                    "gene_symbol")
colnames(luad) <- c("ncbi_gene_id","LUAD_fdr_adjusted_p_value", "LUAD_cancer_sample_med",
                    "LUAD_normal_sample_med", "LUAD_log2_fold_change", "LUAD_p_value", 
                    "gene_symbol")

# merge paad and kirc
paad_kirc <- merge(paad, kirc, by=c("ncbi_gene_id", "gene_symbol"))
# merge paad and kirc and brca
paad_kirc_brca <- merge(paad_kirc, brca, by=c("ncbi_gene_id", "gene_symbol"))
# paad and kirc and brca and luad
all <- merge(paad_kirc_brca, luad, by=c("ncbi_gene_id", "gene_symbol"))

View(all)

# write data to the csv file
#write.csv(all, "all.csv", row.names = FALSE)

# read in data from the csv file
#a <- read.csv("all.csv", header=TRUE)

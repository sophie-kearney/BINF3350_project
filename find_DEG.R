######
# LOAD PACKAGES
######
suppressMessages(library(dplyr))
suppressMessages(library(janitor))

######
# DATA PRE-PROCESSING
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

# remove cancer type column
paad <- paad_p_val[-1]
kirc <- kirc_p_val[-1]
brca <- brca_p_val[-1]
luad <- luad_p_val[-1]


head(brca)
head(kirc)
head(luad)
head(paad)

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

# extracting relevant columns
brcavp <- subset(brca, select = c('ncbi_gene_id', 'BRCA_fdr_adjusted_p_value', 'BRCA_log2_fold_change'))
kircvp <- subset(kirc, select = c('ncbi_gene_id', 'KIRC_fdr_adjusted_p_value', 'KIRC_log2_fold_change'))
luadvp <- subset(luad, select = c('ncbi_gene_id', 'LUAD_fdr_adjusted_p_value', 'LUAD_log2_fold_change'))
paadvp <- subset(paad, select = c('ncbi_gene_id', 'PAAD_fdr_adjusted_p_value', 'PAAD_log2_fold_change'))

View(brcavp)
View(kircvp)
View(luadvp)
View(paadvp)

######
# VOLCANO PLOT
######

suppressMessages(library(ggplot2))

# estatblish cutoff values
# padj.cutoff <- 0.05
# lfc.cutoff <- 0.58

# BRCA Volcano Plot
ggplot(data = brcavp, aes(x = BRCA_log2_fold_change, y = -log10(BRCA_fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("BRCA Volcano Plot")

# KIRC Volcano Plot
ggplot(data = kircvp, aes(x = KIRC_log2_fold_change, y = -log10(KIRC_fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("KIRC Volcano Plot")

# LUAD Volcano Plot
ggplot(data = luadvp, aes(x = LUAD_log2_fold_change, y = -log10(LUAD_fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("LUAD Volcano Plot")

# PAAD Volcano Plot
ggplot(data = paadvp, aes(x = PAAD_log2_fold_change, y = -log10(PAAD_fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("PAAD Volcano Plot")



hist(brcavp$BRCA_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="BRCA - frequencies of adj. p-values\n(all genes)")

hist(kircvp$KIRC_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="KIRC - frequencies of adj. p-values\n(all genes)")

hist(luadvp$LUAD_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="LUAD - frequencies of adj. p-values\n(all genes)")

hist(paadvp$PAAD_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="PAAD - frequencies of adj. p-values\n(all genes)")


######
# SUBSET LUAD
######

luad_deg <- luad[( ((luad$LUAD_log2_fold_change > 1.5) 
                 | (luad$LUAD_log2_fold_change < -1.5))
                 & (-log10(luad$LUAD_fdr_adjusted_p_value) > 10) ),]

test <- luadvp %>% mutate(Color = ifelse((((luad$LUAD_log2_fold_change > 1.5) 
                                          | (luad$LUAD_log2_fold_change < -1.5))
                                         & (-log10(luad$LUAD_fdr_adjusted_p_value) > 10)), 
                                         "black", "red")) 

ggplot(data = test, aes(x = LUAD_log2_fold_change, 
                          y = -log10(LUAD_fdr_adjusted_p_value),
                          color=Color)) + geom_point(size = 0.5) + 
       coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
       ggtitle("LUAD Volcano Plot") + 
       scale_color_manual(values = c("red", "black")) + 
       labs(y= "-log10(FDR Adjusted P Value)", x = "Log Fold Change") + 
       theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


########################################################################
########################################################################
########################################################################



######
# FIND COMMON GENES 
######

# create unique column names for each cancer
# colnames(paad) <- c("ncbi_gene_id","PAAD_fdr_adjusted_p_value", "PAAD_cancer_sample_med",
                    # "PAAD_normal_sample_med", "PAAD_log2_fold_change", "PAAD_p_value", 
                    # "gene_symbol")
# colnames(kirc) <- c("ncbi_gene_id","KIRC_fdr_adjusted_p_value", "KIRC_cancer_sample_med",
                    # "KIRC_normal_sample_med", "KIRC_log2_fold_change", "KIRC_p_value", 
                    # "gene_symbol")
# colnames(brca) <- c("ncbi_gene_id","BRCA_fdr_adjusted_p_value", "BRCA_cancer_sample_med",
                    # "BRCA_normal_sample_med", "BRCA_log2_fold_change", "BRCA_p_value", 
                    # "gene_symbol")
# colnames(luad) <- c("ncbi_gene_id","LUAD_fdr_adjusted_p_value", "LUAD_cancer_sample_med",
                    # "LUAD_normal_sample_med", "LUAD_log2_fold_change", "LUAD_p_value", 
                    # "gene_symbol")

# merge paad and kirc
# paad_kirc <- merge(paad, kirc, by=c("ncbi_gene_id", "gene_symbol"))
# merge paad and kirc and brca
# paad_kirc_brca <- merge(paad_kirc, brca, by=c("ncbi_gene_id", "gene_symbol"))
# paad and kirc and brca and luad
# all <- merge(paad_kirc_brca, luad, by=c("ncbi_gene_id", "gene_symbol"))

# View(all)

# write data to the csv file
#write.csv(all, "all.csv", row.names = FALSE)

# read in data from the csv file
# a <- read.csv("all.csv", header=TRUE)

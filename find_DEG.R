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


# extracting relevant columns
brcavp <- subset(brca, select = c('ncbi_gene_id', 'fdr_adjusted_p_value', 'log2_fold_change'))
kircvp <- subset(kirc, select = c('ncbi_gene_id', 'fdr_adjusted_p_value', 'log2_fold_change'))
luadvp <- subset(luad, select = c('ncbi_gene_id', 'fdr_adjusted_p_value', 'log2_fold_change'))
paadvp <- subset(paad, select = c('ncbi_gene_id', 'fdr_adjusted_p_value', 'log2_fold_change'))

View(brcavp)
View(kircvp)
View(luadvp)
View(paadvp)

######
# VOLCANO PLOT
######

suppressMessages(library(ggplot2))

# establish cutoff values

# BRCA Volcano Plot
ggplot(data = brcavp, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("BRCA Volcano Plot") + geom_vline(xintercept = 1.5, color = "red") + geom_vline(xintercept = -1.5, color = "red") + geom_hline(yintercept = 5, color = "red")

# BRCA Colored Volcano Plot
ggplot(brcavp) +
  geom_point(aes(x=log2_fold_change, y=-log10(fdr_adjusted_p_value), color=threshold)) +
  ggtitle("BRCA Volcano Plot Colored") +
  xlab("log2 fold change") + 
  ylab("-log10 FDR adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

# BRCA DEGs
fdr_adjusted_p_value_cutoff <- 0.1
log2_fold_change_cutoff <- 1

threshold <- brcavp$fdr_adjusted_p_value < fdr_adjusted_p_value_cutoff & abs(brcavp$log2_fold_change) > log2_fold_change_cutoff

  
brca_sig <- filter(brca, (brca$log2_fold_change >= 2 | brca$log2_fold_change <= -2) & (-log10(brca$fdr_adjusted_p_value)))

ggplot(data = brcavp, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-1, 1), ylim = c(10, 20))

# KIRC Volcano Plot
ggplot(data = kircvp, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("KIRC Volcano Plot")

# LUAD Volcano Plot
ggplot(data = luadvp, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("LUAD Volcano Plot")

# PAAD Volcano Plot
ggplot(data = paadvp, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + ggtitle("PAAD Volcano Plot")


######
# DEG HISTOGRAMS
######

hist(brcavp$BRCA_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="BRCA - frequencies of adj. p-values\n(all genes)")

hist(kircvp$KIRC_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="KIRC - frequencies of adj. p-values\n(all genes)")

hist(luadvp$LUAD_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="LUAD - frequencies of adj. p-values\n(all genes)")

hist(paadvp$PAAD_fdr_adjusted_p_value, 
     col="grey", border="white", xlab="", ylab="", main="PAAD - frequencies of adj. p-values\n(all genes)")


########################################################################
########################################################################
########################################################################



######
# FIND COMMON GENES 
######

# create unique column names for each cancer
# colnames(paad) <- c("ncbi_gene_id","PAAD_fdr_adjusted_p_value", "PAAD_cancer_sample_med",
                    #"PAAD_normal_sample_med", "PAAD_log2_fold_change", "PAAD_p_value", 
                    #"gene_symbol")
# colnames(kirc) <- c("ncbi_gene_id","KIRC_fdr_adjusted_p_value", "KIRC_cancer_sample_med",
                    #"KIRC_normal_sample_med", "KIRC_log2_fold_change", "KIRC_p_value", 
                    #"gene_symbol")
# colnames(brca) <- c("ncbi_gene_id","BRCA_fdr_adjusted_p_value", "BRCA_cancer_sample_med",
                    #"BRCA_normal_sample_med", "BRCA_log2_fold_change", "BRCA_p_value", 
                    #"gene_symbol")
# colnames(luad) <- c("ncbi_gene_id","LUAD_fdr_adjusted_p_value", "LUAD_cancer_sample_med",
                    #"LUAD_normal_sample_med", "LUAD_log2_fold_change", "LUAD_p_value", 
                    #"gene_symbol")

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
#a <- read.csv("all.csv", header=TRUE)

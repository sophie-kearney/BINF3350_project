if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

# This is required to make one of the plots that clusterProfiler will make
if (!("ggupset" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ggupset", update = FALSE)
}

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}

if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

brca <- read.csv(file='data/brca_sig.csv')
coad <- read.csv(file='data/coad_sig.csv')
kirc <- read.csv(file='data/kirc_sig.csv')
luad <- read.csv(file='data/luad_sig.csv')
shared <- read.csv(file='data/sharedDEGs.csv')



# reading in input from deseq2
df = read.csv("drosphila_example_de.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- brca$log2_fold_change

# name the vector
names(original_gene_list) <- brca$cancer_type

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = brca

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2_fold_change

# Name the vector
names(genes) <- sig_genes_df$cancer_type

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]

go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'GO',
                      readable = True,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


####
# edits
####

library(clusterProfiler)

brca <- read.csv(file='cross-comparison/brca_diff.csv')
coad <- read.csv(file='cross-comparison/coad_diff.csv')
kirc <- read.csv(file='cross-comparison/kirc_diff.csv')
luad <- read.csv(file='cross-comparison/luad_diff.csv')
shared <- read.csv(file='cross-comparison/shared_DEGs.csv')

luad_ego <- enrichGO(gene     = luad$ncbi_gene_id,
                     OrgDb         = "org.Hs.eg.db",
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

luad_results <-  luad_ego@result
rownames(luad_results) <- 1:nrow(luad_results)
View(luad_results)

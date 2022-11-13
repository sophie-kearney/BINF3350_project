######
# ENVIRONMENT SET-UP
######
suppressMessages(library(dplyr))
suppressMessages(library(janitor))
suppressMessages(library(ggplot2))


######
# DATA PRE-PROCESSING
######
brca_deg <- read.csv(file = 'data/BRCA.csv', sep = '\t')
brca_deg <- brca_deg %>%
  clean_names()
head(brca_deg)

coad_deg <- read.csv(file = 'data/COAD.csv', sep = '\t')
coad_deg <- coad_deg %>%
  clean_names()
head(coad_deg)

kirc_deg <- read.csv(file = 'data/KIRC.csv', sep = '\t')
kirc_deg <- kirc_deg %>%
  clean_names()
head(kirc_deg)

luad_deg <- read.csv(file = 'data/LUAD.csv', sep = '\t')
luad_deg <- luad_deg %>%
  clean_names()
head(luad_deg)


######
# VOLCANO PLOTS (ALL DEGs)
######
ggplot(data = brca_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + ggtitle("BRCA, Volcano Plot")

ggplot(data = coad_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + ggtitle("COAD, Volcano Plot")

ggplot(data = kirc_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + ggtitle("KIRC, Volcano Plot")

ggplot(data = luad_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + geom_point(size = 0.5) + coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + ggtitle("LUAD, Volcano Plot")


######
# EXTRACTING MOST SIGNIFICANT DEGs FOR EACH CANCER TYPE
######

# BRCA
brca_sig <- brca_deg[(((brca_deg$log2_fold_change > 1.5) 
                       | (brca_deg$log2_fold_change < -1.5))
                      & (-log10(brca_deg$fdr_adjusted_p_value) > 10) ),]

write.csv(brca_sig, "brca_sig.csv", row.names = FALSE)

# COAD
coad_sig <- coad_deg[(((coad_deg$log2_fold_change > 1.5) 
                       | (coad_deg$log2_fold_change < -1.5))
                      & (-log10(coad_deg$fdr_adjusted_p_value) > 10) ),]

write.csv(coad_sig, "coad_sig.csv", row.names = FALSE)

# KIRC
kirc_sig <- kirc_deg[(((kirc_deg$log2_fold_change > 1.5) 
                       | (kirc_deg$log2_fold_change < -1.5))
                      & (-log10(kirc_deg$fdr_adjusted_p_value) > 10) ),]

write.csv(kirc_sig, "kirc_sig.csv", row.names = FALSE)

# LUAD
luad_sig <- luad_deg[(((luad_deg$log2_fold_change > 1.5) 
                       | (luad_deg$log2_fold_change < -1.5))
                      & (-log10(luad_deg$fdr_adjusted_p_value) > 10) ),]

write.csv(luad_sig, "luad_sig.csv", row.names = FALSE)


######
# ANNOTATED VOLCANO PLOTS (ALL DEGs, WITH MOST SIGNIFICANT DEGs HIGHLIGHTED)
######
# BRCA
ggplot(data = brca_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  ggtitle("BRCA, Annotated Volcano Plot") + 
  geom_point(pch = 16, size = 0.5) +
  geom_point(data = brca_sig, 
             aes(x = log2_fold_change,y = -log10(fdr_adjusted_p_value)), 
             color = 'red',
             pch = 16,
             size = 0.5)

# COAD
ggplot(data = coad_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  ggtitle("COAD, Annotated Volcano Plot") + 
  geom_point(pch = 16, size = 0.5) +
  geom_point(data = coad_sig, 
             aes(x = log2_fold_change,y = -log10(fdr_adjusted_p_value)), 
             color = 'red',
             pch = 16,
             size = 0.5)

# KIRC
ggplot(data = kirc_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  ggtitle("KIRC, Annotated Volcano Plot") + 
  geom_point(pch = 16, size = 0.5) +
  geom_point(data = kirc_sig, 
             aes(x = log2_fold_change,y = -log10(fdr_adjusted_p_value)), 
             color = 'red',
             pch = 16,
             size = 0.5)

# LUAD
ggplot(data = luad_deg, aes(x = log2_fold_change, y = -log10(fdr_adjusted_p_value))) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  ggtitle("LUAD, Annotated Volcano Plot") + 
  geom_point(pch = 16, size = 0.5) +
  geom_point(data = luad_sig, 
             aes(x = log2_fold_change,y = -log10(fdr_adjusted_p_value)), 
             color = 'red',
             pch = 16,
             size = 0.5)


######
# VOLCANO PLOTS (MOST SIGNIFICANT DEGs ONLY)
######

# BRCA
brca_vp <- brca_sig %>% mutate(Color = ifelse((((brca_sig$log2_fold_change > 1.5) 
                                                | (brca_sig$log2_fold_change < -1.5))
                                               & (-log10(brca_sig$fdr_adjusted_p_value) > 10)), 
                                              "black", "red")) 

ggplot(data = brca_vp, aes(x = log2_fold_change, 
                           y = -log10(fdr_adjusted_p_value),
                           color=Color)) + geom_point(pch = 16, size = 0.5) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  ggtitle("BRCA, Most Significant DEGs") + 
  scale_color_manual(values = c("red", "black")) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# COAD
coad_vp <- coad_sig %>% mutate(Color = ifelse((((coad_sig$log2_fold_change > 1.5) 
                                                | (coad_sig$log2_fold_change < -1.5))
                                               & (-log10(coad_sig$fdr_adjusted_p_value) > 10)), 
                                              "black", "red")) 

ggplot(data = coad_vp, aes(x = log2_fold_change, 
                           y = -log10(fdr_adjusted_p_value),
                           color=Color)) + geom_point(pch = 16, size = 0.5) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  ggtitle("COAD, Most Significant DEGs") + 
  scale_color_manual(values = c("red", "black")) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# KIRC
kirc_vp <- kirc_sig %>% mutate(Color = ifelse((((kirc_sig$log2_fold_change > 1.5) 
                                                | (kirc_sig$log2_fold_change < -1.5))
                                               & (-log10(kirc_sig$fdr_adjusted_p_value) > 10)), 
                                              "black", "red")) 

ggplot(data = kirc_vp, aes(x = log2_fold_change, 
                           y = -log10(fdr_adjusted_p_value),
                           color=Color)) + geom_point(pch = 16, size = 0.5) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  ggtitle("KIRC, Most Significant DEGs") + 
  scale_color_manual(values = c("red", "black")) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# LUAD
luad_vp <- luad_sig %>% mutate(Color = ifelse((((luad_sig$log2_fold_change > 1.5) 
                                                | (luad_sig$log2_fold_change < -1.5))
                                               & (-log10(luad_sig$fdr_adjusted_p_value) > 10)), 
                                              "black", "red")) 

ggplot(data = luad_vp, aes(x = log2_fold_change, 
                           y = -log10(fdr_adjusted_p_value),
                           color=Color)) + geom_point(pch = 16, size = 0.5) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) + 
  ggtitle("LUAD, Most Significant DEGs") + 
  scale_color_manual(values = c("red", "black")) + 
  labs(y= "-log10(FDR-Adjusted P-Value)", x = "Log2-Fold Change") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))












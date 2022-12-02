suppressMessages(library(janitor))
suppressMessages(library(dplyr))
suppressMessages(library(EnhancedVolcano))

brca <- read.table(file='data/brca.csv', header=TRUE, sep='\t')
brca <- brca %>%
  clean_names()

coad <- read.table(file='data/coad.csv', header=TRUE, sep='\t')
coad <- coad %>%
  clean_names()

kirc <- read.table(file='data/kirc.csv', header=TRUE, sep='\t')
kirc <- kirc %>%
  clean_names()

luad <- read.table(file='data/luad.csv', header=TRUE, sep='\t')
luad <- luad %>%
  clean_names()

sigShared = c('IQGAP3', 'RRM2', 'BUB1', 'COL5A2', 'FAM107A', 'ADH1B', 'TRIP13', 
              'SEMA6A', 'KIF20A', 'SOSTDC1', 'ANLN', 'CD36', 'ADAMDEC1', 
              'CTHRC1', 'MELK', 'COL5A1', 'CEP55', 'SCN4B', 'NDRG2', 'SEMA6D', 
              'TNS4', 'COL1A1', 'CA4', 'ABCA8', 'BIRC5', 'PRR36', 'TPX2', 
              'MYBL2', 'UBE2C', 'CDC45', 'MMP11')


#BRCA
keyvals.brca.color <- ifelse (
  brca$log2_fold_change < -1.5 & -log10(brca$fdr_adjusted_p_value) > 10, 'seagreen2',
  ifelse(brca$log2_fold_change > 1.5 & -log10(brca$fdr_adjusted_p_value) > 10, 'maroon2',
         'black'))
keyvals.brca.color[is.na(keyvals.brca.color)] <- 'black'
names(keyvals.brca.color)[keyvals.brca.color == 'maroon2'] <- 'Up-Regulated'
names(keyvals.brca.color)[keyvals.brca.color == 'seagreen2'] <- 'Down-Regulated'
names(keyvals.brca.color)[keyvals.brca.color == 'black'] <- 'N/A'

keyvals.brca.shape <- ifelse (
  brca$gene_symbol %in% sigShared, 12, 19)
keyvals.brca.shape[is.na(keyvals.brca.shape)] <- 19
names(keyvals.brca.shape)[keyvals.brca.shape == 12] <- 'Shared Gene'
names(keyvals.brca.shape)[keyvals.brca.shape == 19] <- 'Non-shared Gene'

EnhancedVolcano(brca,
                lab = brca$gene_symbol,
                x = 'log2_fold_change',
                y = 'fdr_adjusted_p_value',
                selectLab = sigShared,
                labFace = 'plain',
                boxedLabels = TRUE,
                pCutoff = 10e-11,
                FCcutoff = 1.5,
                labSize = 3,
                cutoffLineType = 'twodash',
                colCustom = keyvals.brca.color,
                shapeCustom = keyvals.brca.shape,
                colAlpha = 0.6,
                legendPosition = 'top',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 40,
                title = 'BRCA Volcano Plot',
                subtitle = '') + ggplot2::coord_cartesian(xlim=c(-8, 8)) +
  ggplot2::scale_x_continuous(breaks=seq(-8,8, 1))



#COAD
keyvals.coad.color <- ifelse (
  coad$log2_fold_change < -1.5 & -log10(coad$fdr_adjusted_p_value) > 10, 'seagreen2',
  ifelse(coad$log2_fold_change > 1.5 & -log10(coad$fdr_adjusted_p_value) > 10, 'maroon2',
         'black'))
keyvals.coad.color[is.na(keyvals.coad.color)] <- 'black'
names(keyvals.coad.color)[keyvals.coad.color == 'maroon2'] <- 'Up-Regulated'
names(keyvals.coad.color)[keyvals.coad.color == 'seagreen2'] <- 'Down-Regulated'
names(keyvals.coad.color)[keyvals.coad.color == 'black'] <- 'N/A'

keyvals.coad.shape <- ifelse (
  coad$gene_symbol %in% sigShared, 12, 19)
keyvals.luad.shape[is.na(keyvals.coad.shape)] <- 19
names(keyvals.coad.shape)[keyvals.coad.shape == 12] <- 'Shared Gene'
names(keyvals.coad.shape)[keyvals.coad.shape == 19] <- 'Non-shared Gene'

EnhancedVolcano(coad,
                lab = coad$gene_symbol,
                x = 'log2_fold_change',
                y = 'fdr_adjusted_p_value',
                selectLab = sigShared,
                labFace = 'plain',
                boxedLabels = TRUE,
                pCutoff = 10e-11,
                FCcutoff = 1.5,
                labSize = 3,
                cutoffLineType = 'twodash',
                colCustom = keyvals.coad.color,
                shapeCustom = keyvals.coad.shape,
                colAlpha = 0.6,
                legendPosition = 'top',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 40,
                title = 'COAD Volcano Plot',
                subtitle = '') + ggplot2::coord_cartesian(xlim=c(-9, 7)) +
  ggplot2::scale_x_continuous(breaks=seq(-9,7, 1))


#KIRC
keyvals.kirc.color <- ifelse (
  kirc$log2_fold_change < -1.5 & -log10(kirc$fdr_adjusted_p_value) > 10, 'seagreen2',
  ifelse(kirc$log2_fold_change > 1.5 & -log10(kirc$fdr_adjusted_p_value) > 10, 'maroon2',
         'black'))
keyvals.kirc.color[is.na(keyvals.kirc.color)] <- 'black'
names(keyvals.kirc.color)[keyvals.kirc.color == 'maroon2'] <- 'Up-Regulated'
names(keyvals.kirc.color)[keyvals.kirc.color == 'seagreen2'] <- 'Down-Regulated'
names(keyvals.kirc.color)[keyvals.kirc.color == 'black'] <- 'N/A'

keyvals.kirc.shape <- ifelse (
  kirc$gene_symbol %in% sigShared, 12, 19)
keyvals.kirc.shape[is.na(keyvals.kirc.shape)] <- 19
names(keyvals.kirc.shape)[keyvals.kirc.shape == 12] <- 'Shared Gene'
names(keyvals.kirc.shape)[keyvals.kirc.shape == 19] <- 'Non-shared Gene'

EnhancedVolcano(kirc,
                lab = kirc$gene_symbol,
                x = 'log2_fold_change',
                y = 'fdr_adjusted_p_value',
                selectLab = sigShared,
                labFace = 'plain',
                boxedLabels = TRUE,
                pCutoff = 10e-11,
                FCcutoff = 1.5,
                labSize = 3,
                cutoffLineType = 'twodash',
                colCustom = keyvals.kirc.color,
                shapeCustom = keyvals.kirc.shape,
                colAlpha = 0.6,
                legendPosition = 'top',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 40,
                title = 'KIRC Volcano Plot',
                subtitle = '') + ggplot2::coord_cartesian(xlim=c(-8, 8)) +
  ggplot2::scale_x_continuous(breaks=seq(-8,8, 1))



#LUAD
keyvals.luad.color <- ifelse (
  luad$log2_fold_change < -1.5 & -log10(luad$fdr_adjusted_p_value) > 10, 'seagreen2',
  ifelse(luad$log2_fold_change > 1.5 & -log10(luad$fdr_adjusted_p_value) > 10, 'maroon2',
         'black'))
keyvals.luad.color[is.na(keyvals.luad.color)] <- 'black'
names(keyvals.luad.color)[keyvals.luad.color == 'maroon2'] <- 'Up-Regulated'
names(keyvals.luad.color)[keyvals.luad.color == 'seagreen2'] <- 'Down-Regulated'
names(keyvals.luad.color)[keyvals.luad.color == 'black'] <- 'N/A'

keyvals.luad.shape <- ifelse (
  luad$gene_symbol %in% sigShared, 12, 19)
keyvals.luad.shape[is.na(keyvals.luad.shape)] <- 19
names(keyvals.luad.shape)[keyvals.luad.shape == 12] <- 'Shared Gene'
names(keyvals.luad.shape)[keyvals.luad.shape == 19] <- 'Non-shared Gene'

EnhancedVolcano(luad,
                lab = luad$gene_symbol,
                x = 'log2_fold_change',
                y = 'fdr_adjusted_p_value',
                selectLab = sigShared,
                labFace = 'plain',
                boxedLabels = TRUE,
                pCutoff = 10e-11,
                FCcutoff = 1.5,
                labSize = 3,
                cutoffLineType = 'twodash',
                colCustom = keyvals.luad.color,
                shapeCustom = keyvals.luad.shape,
                colAlpha = 0.6,
                legendPosition = 'top',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 40,
                title = 'LUAD Volcano Plot',
                subtitle = '') + ggplot2::coord_cartesian(xlim=c(-8, 8)) +
                ggplot2::scale_x_continuous(breaks=seq(-8,8, 1))

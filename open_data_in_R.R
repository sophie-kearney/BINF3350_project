library(dplyr)
library(janitor)
lung <- read.table(file="/Users/sophiekearney/Desktop/LUAD_Differential_Gene_Expression_Table.txt", 
                   header=TRUE, sep="\t")
lung <- lung %>%
  clean_names()


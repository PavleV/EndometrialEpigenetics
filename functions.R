# load libraries

library(tidyverse)
library(ggplot2)
library(GenomicRanges)

# load data

RNA_data_df <- readRDS("./Data/ALL_avg_celltype_RNA.rds")
Peak_data_df <- readRDS("./Data/ALL_avg_celltype_peaks.rds")

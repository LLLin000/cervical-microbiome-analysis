library(devtools)
library(BiocManager)
library(tidytree)
library(igraph)
library(ggplot2)
library(magrittr)
library(microeco)

OTU <- read.csv("OTU.csv",header=T, check.names=FALSE, row.names=1)
Sample <- read.csv("Metadata.csv", header=T, check.names=FALSE, row.names=1)
Sample$Sample <- rownames(Sample)
Tax <- read.csv("Tax.csv", header=T, check.names=FALSE, row.names=1)

Tax = microeco::tidy_taxonomy(Tax)

MicrotableObj <- microeco::microtable$new(sample_table = Sample,
                                          otu_table = OTU,
                                          tax_table = Tax,
                                          auto_tidy = F)

LEfSeObj <- microeco::trans_diff$new(dataset = MicrotableObj,
                                     method = "lefse",
                                     group = "SampleType",
                                     lefse_subgroup = NULL,
                                     alpha = 0.05,
                                     p_adjust_method = "none")

write.csv(LEfSeObj$res_diff,"LEfSe_Diff.csv")

LEfSeObj$plot_diff_bar(threshold = 2.0,
                       width = 0.8,
                       group_order = NULL)

LEfSeObj$plot_diff_cladogram(#use_taxa_num = 100,
                             # use_feature_num = 40, 
                             clade_label_level = 4, 
                             group_order = NULL,
                             color = RColorBrewer::brewer.pal(8, "Dark2"),
                             select_show_labels = NULL,
                             only_select_show = FALSE,
                             sep = "|",
                             branch_size = 0.5,
                             alpha = 0.2,
                             clade_label_size = 1.5)

########## tiff ##########
tiff('LEfSe_Cladogram.tiff',res=300,width = 3000,height = 3000,units='px')
LEfSeObj$plot_diff_cladogram(#use_taxa_num = 100,
                             # use_feature_num = 40, 
                             clade_label_level = 4, 
                             group_order = NULL,
                             color = RColorBrewer::brewer.pal(8, "Dark2"),
                             select_show_labels = NULL,
                             only_select_show = FALSE,
                             sep = "|",
                             branch_size = 0.5,
                             alpha = 0.2,
                             clade_label_size = 1.5)
dev.off()

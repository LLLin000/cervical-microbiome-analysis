library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(pheatmap)
library(purrr)
library(dplyr)
library(tidyr)
library(DT)

otu_data <- read.csv("OTU_TaxID.csv", check.names = FALSE)
group_data <- read.csv("CINGroup.csv", row.names = 1)
ref <- "NC"

rownames(otu_data) <- otu_data$`#NAME`
otu_data <- otu_data[, -(1:2)]


common_samples <- intersect(colnames(otu_data), rownames(group_data)) 
otu_data <- otu_data[, common_samples, drop = FALSE]
group_data <- group_data[common_samples, , drop = FALSE]

group_data <- group_data[colnames(otu_data), , drop = FALSE]
group_data$Group <- relevel(factor(group_data$Group), ref = ref)

otu_ps <- otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)
sample_ps <- sample_data(group_data)
ps <- phyloseq(otu_ps, sample_ps)

out <- ancombc(
  data = ps,
  formula = "Group", 
  p_adj_method = "holm", 
  prv_cut = 0.10,
  lib_cut = 1000,  
  group = "Group",
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = TRUE
)

res = out$res
res_global = out$res_global

col_name = colnames(res$lfc)
# ref = sub("^Group", "", setdiff(paste0("Group",group_data$Group|>unique()),colnames(out$res$lfc)[-2:0]))

write.csv(res_global, "./res_global.csv", row.names = FALSE)

tab_lfc = res$lfc
colnames(tab_lfc) = col_name
write.csv(tab_lfc, "./tab_lfc.csv", row.names = FALSE)

tab_se = res$se
colnames(tab_se) = col_name
write.csv(tab_se, "./tab_se.csv", row.names = FALSE)

tab_w = res$W
colnames(tab_w) = col_name
write.csv(tab_w, "./tab_w.csv", row.names = FALSE)

tab_p = res$p_val
colnames(tab_p) = col_name
write.csv(tab_p, "./tab_p.csv", row.names = FALSE)

tab_q = res$q
colnames(tab_q) = col_name
write.csv(tab_q, "./tab_q.csv", row.names = FALSE)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
write.csv(tab_diff, "./tab_diff.csv", row.names = FALSE)

merge_res <- function(res, element) {
  df <- res[[element]]
  if (is.null(df)) return(NULL)
  df_long <- pivot_longer(df, -taxon, names_to = "comparison", values_to = element)
  return(df_long)
}

df_lfc  <- merge_res(res, "lfc")
df_q    <- merge_res(res, "q_val")
df_diff <- merge_res(res, "diff_abn")
res_long <- reduce(list(df_lfc, df_q, df_diff), 
                   left_join, by = c("taxon", "comparison"))

plot_heatmap <- function(df, q_cutoff = NULL, out_file = "./heatmap.svg") {
  if(is.null(q_cutoff))
  {
    df_sig <- df
  }
  else
  {
    df_sig <- df %>% filter(q_val < q_cutoff)
  }
  
  mat <- df_sig %>% 
    filter(comparison != "(Intercept)") %>% 
    select(taxon, comparison, lfc) %>%
    pivot_wider(names_from = comparison, values_from = lfc) %>%
    as.data.frame()
  
  rownames(mat) <- mat$taxon
  mat <- mat[, -1, drop = FALSE]
  
  svg(out_file, width = 16, height = 24)
  pheatmap(mat,
           cluster_rows = TRUE, cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = paste0("Log Fold Changes (Ref: ",ref,")"),
           display_numbers = TRUE,     
           number_format = "%.2f",    
           fontsize_number = 10    
           )
  dev.off()
}

plot_heatmap(res_long, q_cutoff = NULL, out_file = "./heatmap.svg")


################################################################################

sig_taxa <- res_global %>%
  filter(diff_abn == TRUE) %>%
  pull(taxon)

mat <- tab_lfc %>%
  filter(taxon %in% sig_taxa) %>%
  select(-`(Intercept)`) %>%
  mutate(
    taxon = factor(taxon, levels = sig_taxa)
  )

rownames(mat) <- mat$taxon
mat <- mat[, -1, drop = FALSE]

svg("./global_heatmap.svg", width = 16, height = 24)
pheatmap(mat,
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = paste0("Log fold changes for globally significant taxa (Ref: ",ref,")"),
         display_numbers = TRUE,   
         number_format = "%.2f",   
         fontsize_number = 10   
)
dev.off()

library(data.table)
library(fgsea)
library(ggplot2)
library(dplyr)

#==========================prepare fgsea analysis===============================
#use human vs clean gut as an example, can change data set later
# read in DEseq result data
library(readxl)
res <- read_excel("Ag_h_gvsAg_c1_g_deg.xlsx")

# confirm important colum
stopifnot(all(c("gene_id", "log2FoldChange", "pvalue", "padj") %in% names(res)))

# only leave protein coding genes
if ("gene_biotype" %in% names(res)) {
  res_use <- res[res$gene_biotype == "protein_coding", ]
} else {
  res_use <- res
}

# remove NA
res_use <- res_use[!is.na(res_use$gene_id) &
                     !is.na(res_use$log2FoldChange) &
                     !is.na(res_use$pvalue), ]


# 5) build ranking（GSEA）
res_use$rank_score <- sign(res_use$log2FoldChange) * (-log10(res_use$pvalue))

# deal with Infinite (since some pvalue are zero, need to replace pvalue =0 to a minimum finite number)
finite_min <- min(res_use$rank_score[is.finite(res_use$rank_score)], na.rm = TRUE)
res_use$rank_score[!is.finite(res_use$rank_score)] <- finite_min

# 6) remove duplicate, only leave the biggest |rank_score| for the same gene_id
res_use <- res_use[order(res_use$gene_id, -abs(res_use$rank_score)), ]
res_use <- res_use[!duplicated(res_use$gene_id), ]

# 7)generate fgsea "named vector"
ranks <- res_use$rank_score
names(ranks) <- res_use$gene_id
ranks <- sort(ranks, decreasing = TRUE)

#export ranked gene list for source data
rank_table <- res_use %>%
  dplyr::select(gene_id, rank_score, log2FoldChange, pvalue) %>%
  arrange(desc(rank_score))

write.csv(rank_table,
          "Fig3_human_gut_GSEA_ranked_gene_list.csv",
          row.names = FALSE)
#finish exporting

#add functional group info
master <- read_excel("novgenmaster_Anopheles_change_update_12232025.xlsx")

library(dplyr)
master0 <- master %>%
  mutate(FG = coalesce(New_assigned_Function, Functional_Group)) %>%
  dplyr::select(gene_id,FG)


master0$FG <- toupper(master0$FG) #turn to uppe-case

#change format of functional group
master1 <- master0 %>%
  mutate(FG=replace(FG, FG=="D", "DIV")) %>%
  mutate(FG=replace(FG, FG=="I(LRR)", "I/A")) %>%
  mutate(FG=replace(FG, FG=="I(IG)", "I/A")) %>%
  mutate(FG=replace(FG, FG=="I(FBN)", "I/A")) %>%
  mutate(FG=replace(FG, FG=="I", "I/A")) %>%
  mutate(FG=replace(FG, FG=="CS", "C/S")) %>%
  mutate(FG=replace(FG, FG=="CSR", "C/S")) %>%
  mutate(FG=replace(FG, FG=="CSR(OLF)", "C/S")) %>%
  mutate(FG=replace(FG, FG=="DIG", "P/D")) %>%
  mutate(FG=replace(FG, FG=="PROT", "P/D")) %>%
  mutate(FG=replace(FG, FG=="REDOX/STRESS/MIT", "R/S/M")) %>%
  mutate(FG=replace(FG, FG=="RSM", "R/S/M")) %>%
  mutate(FG=replace(FG, FG=="RTT", "R/T/T")) %>%
  mutate(FG=replace(FG, FG=="M", "MET")) %>%
  mutate(FG=replace(FG, FG=="U", "UKN"))

#check if there is extra Functional group that I did not change format, should be 9
n_distinct(master1$FG) 

#build gene set for GSEA
gene_sets <- master1 %>%
  filter(!is.na(FG)) %>%
  group_by(FG) %>%
  summarise(genes = list(unique(gene_id))) %>%
  tibble::deframe()

#check the size of genes_sets for proper number of genes
sort(lengths(gene_sets))

#check if rank vector and genes_sets name are the aligned
sum(res_use$gene_id %in% unique(unlist(gene_sets)))
sum(res_use$gene_id %in% unique(unlist(gene_sets))) / length(unique(unlist(gene_sets)))
length(intersect(names(ranks), unique(unlist(gene_sets))))

# check if one gene_id has multiple functions
tab <- table(unlist(gene_sets))
sum(tab > 1)

#=======================================above are preparation work for fgsea=================================
# only keep the gene_id in ranks
gene_sets2 <- lapply(gene_sets, function(x) intersect(x, names(ranks)))
gene_sets2 <- gene_sets2[lengths(gene_sets2) >= 15]  # filter minSize=15

fg <- fgseaMultilevel(pathways = gene_sets2, stats = ranks,
                      minSize = 15, maxSize = 2000) #minSize and maxSize depend on the size of functional group, leave DIV and UNK out
fg <- fg %>% arrange(padj)
fg

# visualization
fg_plot <- fg %>%
  mutate(pathway = factor(pathway, levels = pathway[order(NES)]))

#plot bar plot
ggplot(fg_plot, aes(x = NES, y = pathway)) +
  geom_col() +
  theme_classic() +
  labs(x = "NES", y = "Functional group")

#plot dot plot
p <- ggplot(fg, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, shape = padj < 0.05)) +
  theme_classic() +
  labs(x = "NES", y = "Functional group", size = "Genes in set")

ggsave(
  filename = "Fig_GSEA_gut_h_c1.eps",
  plot = p,
  device = cairo_ps,
  width = 6,
  height = 5,
  units = "in"
)


#make enrichment curve for each functional group and put them into the supplementary figure
p_CS   <- plotEnrichment(gene_sets2[["C/S"]],  ranks) + ggtitle("C/S")
p_MET  <- plotEnrichment(gene_sets2[["MET"]],  ranks) + ggtitle("MET")
p_RSM  <- plotEnrichment(gene_sets2[["R/S/M"]],ranks) + ggtitle("R/S/M")
p_PD   <- plotEnrichment(gene_sets2[["P/D"]],  ranks) + ggtitle("P/D")
p_UKN  <- plotEnrichment(gene_sets2[["UKN"]],  ranks) + ggtitle("UKN")
p_RTT  <- plotEnrichment(gene_sets2[["R/T/T"]],ranks) + ggtitle("R/T/T")
p_TRP  <- plotEnrichment(gene_sets2[["TRP"]],  ranks) + ggtitle("TRP")
p_IA   <- plotEnrichment(gene_sets2[["I/A"]],  ranks) + ggtitle("I/A")
p_DIV  <- plotEnrichment(gene_sets2[["DIV"]],  ranks) + ggtitle("DIV")

library(patchwork)

p_all <- (p_PD | p_RSM | p_MET) /
  (p_RTT | p_IA | p_TRP) /
  (p_CS | p_DIV  | p_UKN)


ggsave(
  "Supplementary_Fig_Sx_GSEA_enrichment_curves.eps",
  plot = p_all,
  device = cairo_ps,
  width = 180, height = 225, units = "mm",
)

ggsave(
  "Supplementary_Fig_Sx_GSEA_enrichment_curves.tiff",
  plot = p_all,
  width = 180, height = 225, units = "mm",
  dpi = 300
)
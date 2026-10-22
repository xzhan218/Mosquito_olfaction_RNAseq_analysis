#GO and KEGG for DEG

library(readxl)
library(org.Ag.eg.db) #anopheles gambiae 
library(AnnotationDbi)
keytypes(org.Ag.eg.db)
library(ggplot2)
library(gprofiler2) #may be more updated
library(dplyr)

#====use the belowing code==================================================================================================================
#====try goprofiler2========================================================================================================================
#====human gut==============
h_g <- read_excel("Ag_h_gvsAg_c1_g_deg_122724.xlsx")

deg <- subset(h_g, !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange))
#2315

deg_up   <- unique(deg$gene_id[deg$log2FoldChange > 0])
deg_down <- unique(deg$gene_id[deg$log2FoldChange < 0])

#====前期filter跟wcgna一致======
# background（最保守的做法：用该比较里“被检测过”的全部基因）
countdata <- read.csv("gene_count_01052025.csv",row.names = "gene_id",header = T)
dim(countdata)
#13845 24

# Pre-filter data
#filter genes with counts less than 10 in more than 90% of samples
countdata$low = apply(countdata[,1:24],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each group

count2<-countdata[-which(countdata$low>22),] #22 is 90% of 30 samples - get rid of count less than 10 in more than 90% of samples

nrow(count2) #10636

#====第二次filter only count comparison samples=============
samples_bg <- c("Ag_hg1","Ag_hg2","Ag_hg5","Ag_hg6",
                "Ag_c1g1","Ag_c1g2","Ag_c1g4","Ag_c1g6")

cts_bg <- count2[, samples_bg, drop = FALSE]

#做一个轻量表达过滤：至少2个样本 count>=10
bg_genes <- rownames(cts_bg)[ rowSums(cts_bg >= 10) >= 2 ]
length(bg_genes)
#9014

deg_up   <- unique(trimws(deg_up))
deg_down <- unique(trimws(deg_down))
bg_genes <- unique(trimws(bg_genes))

org <- "agambiae"

gost_up <- gost(
  query = deg_up,
  organism = org,
  custom_bg = bg_genes,
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  correction_method = "fdr", #p value adjustment methods
  sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC") #source included BP.MF.CC, KEGG
)

gost_down <- gost(
  query = deg_down,
  organism = org,
  custom_bg = bg_genes,
  significant = TRUE,
  correction_method = "fdr",
  sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC")
)

up_df   <- gost_up$result
down_df <- if (is.null(gost_down)) NULL else gost_down$result

head(up_df[, c("source","term_name","p_value","intersection_size","term_size")])
head(down_df[, c("source","term_name","p_value","intersection_size","term_size")])

up_df[] <- lapply(up_df, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})

down_df[] <- lapply(down_df, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})

write.csv(up_df, "GO_human_gut_up_padj.csv", row.names = FALSE)
write.csv(down_df, "GO_human_gut_down_padj.csv", row.names = FALSE)

#GO
#up plot
plot_df <- up_df %>%
  filter(source == "GO:BP") %>%
  filter(!term_name %in% c(
    "biological_process",
    "cellular process",
    "metabolic process"
  )) %>%
  mutate(
    minus_log10_FDR = -log10(p_value),
    gene_ratio = intersection_size / term_size
  ) %>%
  arrange(p_value)

plot_df_top <- plot_df %>% slice(1:12) #top 12 BP; can change to 8 or 10

p <- ggplot(plot_df_top,
            aes(x = intersection_size,
                y = reorder(term_name, intersection_size))) +
  geom_point(aes(size = gene_ratio,
                 color = minus_log10_FDR)) +
  scale_color_viridis_c(name = "-log10(FDR)") +
  scale_size(range = c(3, 8), name = "Gene ratio") +
  labs(
    x = "Number of genes",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

p <- p +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.25))) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(6, 20, 6, 10)) +
  guides(
    color = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  )

ggsave(
  "GO_up_midgut_human_vs_clean_padj.png",
  p,
  width = 6.8,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "GO_up_midgut_human_vs_clean.eps",
  p,
  device = cairo_ps,
  width = 6.8,
  height = 4.5,
  units = "in",
)


#down plot
plot_df <- down_df %>%
  filter(source == "GO:BP") %>%
  filter(!term_name %in% c(
    "biological_process",
    "cellular process",
    "metabolic process"
  )) %>%
  mutate(
    minus_log10_FDR = -log10(p_value),
    gene_ratio = intersection_size / term_size
  ) %>%
  arrange(p_value)

plot_df_top <- plot_df %>% slice(1:12) #top 12 BP; can change to 8 or 10

p <- ggplot(plot_df_top,
            aes(x = intersection_size,
                y = reorder(term_name, intersection_size))) +
  geom_point(aes(size = gene_ratio,
                 color = minus_log10_FDR)) +
  scale_color_viridis_c(name = "-log10(FDR)") +
  scale_size(range = c(3, 8), name = "Gene ratio") +
  labs(
    x = "Number of genes",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

p <- p +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.25))) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(6, 20, 6, 10)) +
  guides(
    color = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  )

ggsave(
  "GO_down_midgut_human_vs_clean_padj.png",
  p,
  width = 6.8,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "GO_down_midgut_human_vs_clean.eps",
  p,
  device = cairo_ps,
  width = 6.8,
  height = 4.5,
  units = "in"
)

#KEGG
#up plot
plot_df <- up_df %>%
  filter(source == "KEGG") %>%
  filter(!term_name %in% c("KEGG root term", "Metabolic pathways"), term_size <= 300) %>% #remove not wanted kegg root term and the termsize too big (一般来说层级太高都显著）)
  mutate(
    minus_log10_FDR = -log10(p_value),
    gene_ratio = intersection_size / term_size
  ) %>%
  arrange(p_value)

plot_df_top <- plot_df %>% slice(1:12) #top 12 BP; can change to 8 or 10

p <- ggplot(plot_df_top,
            aes(x = intersection_size,
                y = reorder(term_name, intersection_size))) +
  geom_point(aes(size = gene_ratio,
                 color = minus_log10_FDR)) +
  scale_color_viridis_c(name = "-log10(FDR)") +
  scale_size(range = c(3, 8), name = "Gene ratio") +
  labs(
    x = "Number of genes",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

p <- p +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.25))) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(6, 20, 6, 10)) +
  guides(
    color = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  )

ggsave(
  "KEGG_up_midgut_human_vs_clean_padj.png",
  p,
  width = 6.8,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "KEGG_up_midgut_human_vs_clean.eps",
  p,
  device = cairo_ps,
  width = 6.8,
  height = 4.5,
  units = "in"
)

#down plot
plot_df <- down_df %>%
  filter(source == "KEGG") %>%
  filter(!term_name %in% c("KEGG root term", "Metabolic pathways"), term_size <= 300) %>% #remove not wanted kegg root term
  mutate(
    minus_log10_FDR = -log10(p_value),
    gene_ratio = intersection_size / term_size
  ) %>%
  arrange(p_value)

plot_df_top <- plot_df %>% slice(1:12) #top 12 BP; can change to 8 or 10

p <- ggplot(plot_df_top,
            aes(x = intersection_size,
                y = reorder(term_name, intersection_size))) +
  geom_point(aes(size = gene_ratio,
                 color = minus_log10_FDR)) +
  scale_color_viridis_c(name = "-log10(FDR)") +
  scale_size(range = c(3, 8), name = "Gene ratio") +
  labs(
    x = "Number of genes",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

p <- p +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.25))) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(6, 20, 6, 10)) +
  guides(
    color = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  )

ggsave(
  "KEGG_down_midgut_human_vs_clean_padj.png",
  p,
  width = 6.8,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "KEGG_down_midgut_human_vs_clean.eps",
  p,
  device = cairo_ps,
  width = 6.8,
  height = 4.5,
  units = "in"
)
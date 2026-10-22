#load packages
library(WGCNA) 
library(dplyr)
library(tidyverse)
library(magrittr)
library(gridExtra)
library(DESeq2)

#------------------------------------------------------------------------------------------------------------
### A. Get data, and data pre-cleaning steps=====================================
#Read in data and create dds data set
countdata <- read.csv("gene_count_01052025.csv",row.names = "gene_id",header = T)
dim(countdata)
coldata <- read.csv("group_01052025.csv",header=TRUE)

# Pre-filter data
#filter genes with counts less than 10 in more than 90% of samples
countdata$low = apply(countdata[,1:24],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each group

count2<-countdata[-which(countdata$low>22),] #22 is 90% of 30 samples - get rid of count less than 10 in more than 90% of samples

nrow(count2) #10636

#build DESeq dataset for data normalization
dds<-DESeqDataSetFromMatrix(countData=count2[,1:24],colData=coldata,design=~1) #since there was a column "low" in the former step, need to remove it here; design = ~1 means not specifying models

#filter data as DESeq analysis, need to be consistant with DESeq analysis; using 30 and 60 here does not make much difference
keep<- rowSums(counts(dds))>=30
dds<-dds[keep,]
nrow(dds) #10636

#Normalize the data using vst()
vsd<-vst(dds,blind=TRUE)
dat=as.data.frame(assay(vsd)) #assay means get the normalized counts, since vsd is a large DESeqDataSet, not only normalized counts, but other info
colnames(dat)<-names(countdata[1:24])


boxplot(dat) #how do expression plots look overall?
#keep all genes
datExpr0 <- as.data.frame(t(dat[1:10636,1:24]))

#Quality control; find gene outliers and sample outliers; this part is conducted after normalization becuase in author tutorial they used already normalied data
gsg <- goodSamplesGenes(datExpr0)
summary(gsg)
gsg$allOK #FALSE meaning that there are outliers

#detect outlier genes
table(gsg$goodGenes) #10636 ALL true
table(gsg$goodSamples) #all 24 samples passed test

#Save expression and trait data for use in the next steps. Note: to prevent covering the old data, change the name of exported files soon
save(datExpr0,coldata,file = "Anopheles_networkConstruction_1.RData")

### B. Network Construction=====================================
# WGCNA try a range of soft thresholds and create a diagnostic plot for us to pick a threshold value (if correlation is below threshold, remove genes)
# allowWGCNAThreads()          # allow multi-threading (optional); Allowing multi-threading with up to 4 threads.

#load data, this commend is stated as important ,cannot omit
options(stringAsFactors = F)
#load the expression and trait data saved before, not I have changed RData name before, so here the name is a bit different
lnames = load(file = "Anopheles_networkConstruction_1change.RData")
lnames

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(
  datExpr0,             # Input data
  powerVector = power,
  networkType = "signed",
  verbose = 5
)

#visualize to pick power
sft.data <- sft$fitIndices #we need to find the best power with highest R.sq and lowest mean.k

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2) #choosing 12 as power

#covert matrix to numeric
datExpr0[] <- sapply(datExpr0, as.numeric)

soft_power <- 12 #define soft power
temp_cor <- cor # force it to use GCNA cor function (fix a namespace conflict issue)
cor <- WGCNA::cor

#Memory estimate w.r.t blocksize
#maxBlockSize default 5000, but it is better to run all genes in one block, my macbook has memory of 16G, so I can have up to 20000. I will set 18000 because this is larger than my gene number
bwnet <- blockwiseModules(datExpr0,
                          maxBlockSize = 18000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix, if positive and negative are important in addition to relatedness, should use "signed"
                          power = soft_power, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
                          mergeCutHeight = 0.25,
                          saveTOMs = TRUE,
                          verbose = 3
)
#took around 1min

#After run
cor <- temp_cor     # Return cor function to original namespace

#Take a look at the modules
# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#Save module assignment and module eigengene info for subsequent analysis
moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs
geneTree = bwnet$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Anopheles_networkConstruction_2.RData")

### C. Related modules to traits================================================
#load data, this commend is stated as important ,cannot omit
options(stringAsFactors = F)
#load the expression and trait data saved before, not I have changed RData name before, so here the name is a bit different
lnames = load(file = "Anopheles_networkConstruction_1change.RData")
lnames
#load the network data saved before
lnames = load(file = "Anopheles_networkConstruction_2change.RData")
lnames

###Create trait file and binarize categorical variables, if the trait is numeric, no need to change
class(coldata)

tissue_bi <- coldata %>%
  mutate(tissue = ifelse(grepl("midgut",tissue),1,0))  %>%
  dplyr::select(3) #select the sixth column

#Originally used this methods to create binary data, but after ranking, control is not showing in the figure, therefore used the next method for binarization
#coldata$treat <- factor(coldata$treat, levels = c("control","human","air"))
#treat_bi <- binarizeCategoricalColumns(coldata$treat,
#                                        includePairwise = F,
#                                        includeLevelVsAll = T,
#                                        minCount =1
#)

#try One-Hot Encoding (Most Common for Nominal Treatments) for binarizing treat data
treat_bi <- coldata %>%
  mutate(
    control = ifelse(treat == "control", 1, 0),
    human = ifelse(treat == "human", 1, 0),
    air = ifelse(treat == "air", 1, 0)
  ) %>%
  select(6,7,8)

traits <- cbind(tissue_bi,treat_bi) #only include binary columns

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes #note here used moduleColors
MEs = orderMEs(MEs0)
#cor(MEs, traits) calculates the correlation coefficient for each pair of columns from MEs and traits
moduleTraitCor = cor(MEs, traits, use = "p"); #"p" (pairwise complete observations): Uses all pairs of samples with non-missing values for the correlation calculation
#calculate the p-values associated with the correlation coefficients in moduleTraitCor. This allows you to assess the statistical significance of the correlations between module eigengenes (MEs) and traits
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#visualization
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot: module-trait relationships
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

####export GS and MM for WGCNA grouping==========================================
#Gene relationship to trait and important modules: gene significance GS and module membership MM
# Define variable weight containing the tissue column of traits; can change tissue to control, human, or air for different tables
air = as.data.frame(traits$air);
names(air) = "air"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, air, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(air), sep="");
names(GSPvalue) = paste("p.GS.", names(air), sep="")

# Create the starting data frame
geneInfo0 = data.frame(
  moduleColor = moduleColors,
  geneTraitSignificance)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, air, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.air));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo_air.csv") #remember to change name

#visualization and follow-up analysis

library(WGCNA) 
library(ggplot2)
library(dplyr)
library(readr)

#load data, this commend is stated as important ,cannot omit
options(stringAsFactors = F)
#load the expression and trait data saved before
lnames = load(file = "Anopheles_networkConstruction_1change.RData")
lnames
#load the network data saved before
lnames = load(file = "Anopheles_networkConstruction_2change.RData")
lnames

#read in geneInfo_control data.frame, to make dot plot (x_axis: Module membership, y_axis: gene significance)
control <- read_csv("geneInfo_control_change.csv")
colnames(control)


#change coldata treat name and level and make a new condition colum
coldata$Condition<- factor(
  coldata$treat,
  levels = c("control", "human", "air"),
  labels = c("Clean air", "Odor stimulation", "Passive exposure")
)

#================================================================
#heatmap
#try One-Hot Encoding (Most Common for Nominal Treatments) for binarizing treat data

traits <- data.frame(
  tissue  = ifelse(grepl("midgut", coldata$tissue), 1, 0),
  clean   = ifelse(coldata$Condition == "Clean air", 1, 0),
  human   = ifelse(coldata$Condition == "Odor stimulation", 1, 0),
  normal  = ifelse(coldata$Condition == "Passive exposure", 1, 0)
)

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes #note here used moduleColors
MEs = orderMEs(MEs0)
#cor(MEs, traits) calculates the correlation coefficient for each pair of columns from MEs and traits
moduleTraitCor = cor(MEs, traits, use = "p"); #"p" (pairwise complete observations): Uses all pairs of samples with non-missing values for the correlation calculation
#calculate the p-values associated with the correlation coefficients in moduleTraitCor. This allows you to assess the statistical significance of the correlations between module eigengenes (MEs) and traits
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#visualization
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot: module-trait relationships
traitLabels <- c(
  "Tissue\n(midgut vs. carcass)",
  "Clean air",
  "Odor stimulation",
  "Passive exposure"
)

#png("Fig_WGCNA_ModuleTraitHeatmap.png", width = 3000, height = 1800, res = 300)

cairo_ps("Fig_WGCNA_Moduletrait_heatmap.eps",
         width = 10,
         height = 6)

par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = traitLabels,
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  zlim = c(-1,1),
  main = "Module–trait relationships"
)
dev.off()
# remember to export

#export cor and p value for source data
cor_df <- as.data.frame(moduleTraitCor)
p_df <- as.data.frame(moduleTraitPvalue)

cor_df$module <- rownames(cor_df)
p_df$module <- rownames(p_df)

source_data <- merge(cor_df, p_df, by="module", suffixes=c("_cor","_pvalue"))

write.csv(source_data,
          "Fig3_WGCNA_module_trait_relationships.csv",
          row.names = FALSE)
#finish exporting

#===========================================================
#for specific module

greenGenes <- control[control$moduleColor == "green",]

colnames(greenGenes)[1] <- "gene_id"

x <- abs(greenGenes$MM.green)
y <- abs(greenGenes$GS.control)

plot(
  x, y,
  xlab = "Module Membership in green module",
  ylab = "Gene Significance for control",
  main = "Module membership vs. gene significance\n(green module)",
  pch = 16,
  col = "darkgreen"
)

correlation <- cor(x, y, use = "p")

legend(
  "topleft",
  legend = paste("cor =", round(correlation, 2)),
  bty = "n"
)

head(moduleColors)
table(moduleColors)[1:10]
colnames(MEs) # I labeled MEs with moduleLabels

#change back to module color
MEs0 <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)
colnames(MEs)

boxplot(MEs$MEgreen ~ coldata$Condition,
        xlab = "Condition",
        ylab = "Module eigengene (green)",
        main = "MEgreen across conditions")

table(coldata$treat, useNA = "ifany") #number of replicates

#stats：see if conditions have significant difference
kruskal.test(MEs$MEgreen ~ coldata$treat)
#Kruskal-Wallis rank sum test

#data:  MEs$MEgreen by coldata$treat
#Kruskal-Wallis chi-squared = 17.565, df = 2, p-value = 0.0001534

#:which comparisons are significant
pairwise.wilcox.test(MEs$MEgreen, coldata$treat, p.adjust.method = "BH")
#Pairwise comparisons using Wilcoxon rank sum exact test 

#data:  MEs$MEgreen and coldata$treat 

#           Clean air Human odor
#Human odor 0.00023   -         
#Normal air 0.00023   0.02813   

#P value adjustment method: BH 

#
fig3d_data <- data.frame(
  sample = rownames(MEs),
  Condition = coldata$Condition,
  tissue = coldata$tissue,
  MEgreen = MEs$MEgreen
)

write.csv(fig3d_data,
          "Fig3d_MEgreen_across_conditions.csv",
          row.names = FALSE)

kw_green <- kruskal.test(MEs$MEgreen ~ coldata$treat)

pw_green <- pairwise.wilcox.test(
  MEs$MEgreen,
  coldata$treat,
  p.adjust.method = "BH"
)

pw_green_df <- as.data.frame(as.table(pw_green$p.value))
colnames(pw_green_df) <- c("group1", "group2", "adjusted_pvalue")

write.csv(pw_green_df,
          "Fig3d_MEgreen_pairwise_stats.csv",
          row.names = FALSE)

kw_green_df <- data.frame(
  test = "Kruskal-Wallis",
  statistic = as.numeric(kw_green$statistic),
  df = as.numeric(kw_green$parameter),
  pvalue = kw_green$p.value
)

write.csv(kw_green_df,
          "Fig3d_MEgreen_overall_stats.csv",
          row.names = FALSE)

#================================================
#remake Green module boxplot and add statistics:
#1. boxplot with extra space
png("Fig_WGCNA_MEgreen_Boxplot.png", width = 1600, height = 1600, res = 300)

cairo_ps("Fig_WGCNA_MEgreen_Boxplot.eps",
         width = 6,
         height = 6)
boxplot(
  MEs$MEgreen ~ coldata$Condition,
  ylim = c(min(MEs$MEgreen), max(MEs$MEgreen) + 0.4),
  xlab = "Condition",
  ylab = "Module eigengene",
  main = "MEgreen across conditions"
)

# 2. y max
y_max <- max(MEs$MEgreen, na.rm = TRUE)

# 3. significance bars
segments(1, y_max + 0.05, 2, y_max + 0.05)
text(1.5, y_max + 0.07, "***")

segments(1, y_max + 0.12, 3, y_max + 0.12)
text(2, y_max + 0.14, "***")

segments(2, y_max + 0.19, 3, y_max + 0.19)
text(2.5, y_max + 0.21, "*")

dev.off()
#===================================================
#check blue module and yellow module
boxplot(
  MEs$MEblue ~ coldata$treat,
  ylim = c(min(MEs$MEblue), max(MEs$MEblue) + 0.4),
  xlab = "Condition",
  ylab = "Module eigengene (blue)",
  main = "MEblue across conditions"
)
kruskal.test(MEs$MEblue ~ coldata$treat)
pairwise.wilcox.test(MEs$MEblue, coldata$treat, p.adjust.method = "BH")

boxplot(MEs$MEblue ~ coldata$tissue,
        xlab="Tissue", ylab="Module eigengene (blue)",
        main="MEblue across tissues")
kruskal.test(MEs$MEblue ~ coldata$tissue)

boxplot(
  MEs$MEyellow ~ coldata$treat,
  ylim = c(min(MEs$MEyellow), max(MEs$MEyellow) + 0.4),
  xlab = "Condition",
  ylab = "Module eigengene (yellow)",
  main = "MEblue across conditions"
)
kruskal.test(MEs$MEyellow ~ coldata$treat)
pairwise.wilcox.test(MEs$MEyellow, coldata$treat, p.adjust.method = "BH")

boxplot(MEs$MEyellow ~ coldata$tissue,
        xlab="Tissue", ylab="Module eigengene (yellow)",
        main="MEyellow across tissues")
kruskal.test(MEs$MEyellow ~ coldata$tissue)
#Kruskal-Wallis rank sum test

#data:  MEs$MEyellow by coldata$tissue
#Kruskal-Wallis chi-squared = 11.213, df = 1, p-value = 0.0008121

#========================================================================
#continue on green module； Check functional group
green <- subset(greenGenes, moduleColor == "green")

colnames(green)[1] <- "gene_id"

anno <- read_excel("gene_function.xlsx")

green_anno <- merge(green, anno, by = "gene_id", all.x = TRUE)

green_anno <- green_anno%>%
  select(gene_id,moduleColor,GS.control,MM.green,p.MM.green,FG,gene_name)

green_pos <- subset(green_anno, moduleColor == "green" & MM.green > 0)
green_neg <- subset(green_anno, moduleColor == "green" & MM.green < 0)

table(green_pos$FG)
table(green_neg$FG)

#========================================================================
#GO enrichment analysis
library(gprofiler2)

#get yellow module gene id
green_ids <- unique(greenGenes$gene_id)  

pos_ids <- unique(green_pos$gene_id)
neg_ids <- unique(green_neg$gene_id)

bg_ids <- colnames(datExpr0) #background id, used datExpr0 (WGCNA gene dataset)

#check BP enrichment first
ego_pos <- gost(
  query = pos_ids,                #set customized query for GO enrichment
  organism = "agambiae",
  sources = c("GO:BP"),
  correction_method = "fdr",
  custom_bg = bg_ids
)

#check BP enrichment first
ego_neg <- gost(
  query = neg_ids,                #set customized query for GO enrichment
  organism = "agambiae",
  sources = c("GO:BP"),
  correction_method = "fdr",
  custom_bg = bg_ids
)

head(ego_pos$result[, c("term_name","p_value","intersection_size")], 10)
head(ego_neg$result[, c("term_name","p_value","intersection_size")], 10)

#export GO enrichment result
pos_out <- ego_pos$result
neg_out <- ego_neg$result

# 把所有 list 列（比如 intersection）转成用逗号连接的字符串
pos_out[] <- lapply(pos_out, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})
neg_out[] <- lapply(neg_out, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})

write.csv(pos_out, "GO_enrichment_green_pos_BP.csv", row.names = FALSE)
write.csv(neg_out, "GO_enrichment_green_neg_BP.csv", row.names = FALSE)

#MF and CC in postive group
ego_pos_MF <- gost(
  query = pos_ids,
  organism = "agambiae",
  sources = c("GO:MF"),
  correction_method = "fdr",
  custom_bg = bg_ids
)

ego_pos_CC <- gost(
  query = pos_ids,
  organism = "agambiae",
  sources = c("GO:CC"),
  correction_method = "fdr",
  custom_bg = bg_ids
)

head(ego_pos_MF$result[, c("term_name","p_value","intersection_size")], 5)
head(ego_pos_CC$result[, c("term_name","p_value","intersection_size")], 5)

pos_MF_out <- ego_pos_MF$result
pos_MF_out[] <- lapply(pos_MF_out, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})
write.csv(pos_MF_out, "GO_enrichment_green_pos_MF.csv", row.names = FALSE)

#===============================================================================
#===============================================================================
#since tissue type is a strong effect in wgcna, we would like to only look at midgut sample and recalculate ME
# 1) subset midgut samples
keep <- grepl("midgut", coldata$tissue)
dat_mid <- datExpr0[keep, ]
col_mid <- coldata[keep, ]

# 2) traits（只在 midgut 内，Condition one-hot 或 factor 都行）
col_mid$Condition <- factor(
  col_mid$treat,
  levels = c("control", "human", "air"),
  labels = c("Clean air", "Odor stimulation", "Passive exposure")
)

# 3) 重新算 MEs（用同一套 moduleColors）
MEs_mid <- moduleEigengenes(dat_mid, moduleColors)$eigengenes
MEs_mid <- orderMEs(MEs_mid)

# 4) 看 yellow
boxplot(MEs_mid$MEyellow ~ col_mid$Condition)

kruskal.test(MEs_mid$MEyellow ~ col_mid$Condition)
#Kruskal-Wallis rank sum test

#data:  MEs_mid$MEyellow by col_mid$Condition
#Kruskal-Wallis chi-squared = 7.2692, df = 2, p-value = 0.02639

pairwise.wilcox.test(MEs_mid$MEyellow, col_mid$Condition, p.adjust.method="BH")
#Pairwise comparisons using Wilcoxon rank sum exact test 

#data:  MEs_mid$MEyellow and col_mid$Condition 

#             Clean air Human odor
#Human odor   0.086     -         
#Normal air   0.086     0.343     

#P value adjustment method: BH 

table(col_mid$Condition) #yellow module figure replicate number

#export source data
fig3e_data <- data.frame(
  sample = rownames(MEs_mid),
  Condition = col_mid$Condition,
  tissue = col_mid$tissue,
  MEyellow = MEs_mid$MEyellow
)

write.csv(
  fig3e_data,
  "Fig3e_MEyellow_midgut_only.csv",
  row.names = FALSE
)
#finish export
#-------------------------------------------------------------------------------

#置换检验: permutation test
library(lmPerm)
fitp <- lmp(ME ~ tissue * condition, data = df)
anova(fitp)
#Analysis of Variance Table

#Response: ME
#                   Df R Sum Sq R Mean Sq Iter Pr(Prob)   
#  tissue            1  0.06634  0.066342 5000   0.0020 **
#  condition         2  0.16442  0.082208 5000   0.0020 **
#  tissue:condition  2  0.05270  0.026352 5000   0.0118 * 
#  Residuals        18  0.71654  0.039808                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Pr: the probability of observing an effect at least as extreme as the one observed, under random permutation of the labels.

#-------------------------------------------------------------------------------
#refine yellow midgut boxplot to green
MEy <- MEs_mid$MEyellow

#png("Fig_WGCNA_MEyellow_gut_Boxplot.png", width = 1600, height = 1600, res = 300)

cairo_ps("Fig_WGCNA_MEyellow_gut_Boxplot.eps",
         width = 6,
         height = 6)

boxplot(
  MEy ~ col_mid$Condition,
  col = "grey85",
  border = "black",
  lwd = 1.2,
  ylim = c(min(MEy), max(MEy) + 0.35),
  xlab = "Condition",
  ylab = "Module eigengene",
  main = "MEyellow across conditions"
)
dev.off()

#===============================================================================
#function and GO enrichment analysis for yellow module
yellowGenes <- control[control$moduleColor == "yellow",]

colnames(yellowGenes)[1] <- "gene_id"

yellow_anno <- merge(yellowGenes, anno, by = "gene_id", all.x = TRUE)

yellow_anno <- yellow_anno%>%
  select(gene_id,moduleColor,GS.control,MM.yellow,p.MM.yellow,FG,gene_name)

yellow_pos <- subset(yellow_anno, moduleColor == "yellow" & MM.yellow > 0)
yellow_neg <- subset(yellow_anno, moduleColor == "yellow" & MM.yellow < 0)
#yellow Go enrichment analysis
#GO enrichment analysis

#get yellow module gene id
yellow_ids <- unique(yellowGenes$gene_id)

#set up background gene id
bg_ids <- colnames(datExpr0)

#GO enrichment analysis
library(gprofiler2)

ego_yellow_BP <- gost(
  query = yellow_ids,
  organism = "agambiae",
  sources = "GO:BP",
  correction_method = "fdr",
  custom_bg = bg_ids
)

head(ego_yellow_BP$result[, c("term_name","p_value","intersection_size")], 10)

pos_y_ids <- unique(yellow_pos$gene_id)
neg_y_ids <- unique(yellow_neg$gene_id)

#check pos enrichment
ego_y_pos <- gost(
  query = pos_y_ids,                #set customized query for GO enrichment
  organism = "agambiae",
  sources = c("GO:BP"),
  correction_method = "fdr",
  custom_bg = bg_ids
)

#check neg enrichment
ego_y_neg <- gost(
  query = neg_y_ids,                #set customized query for GO enrichment
  organism = "agambiae",
  sources = c("GO:BP"),
  correction_method = "fdr",
  custom_bg = bg_ids
)

#export GO enrichment result
pos_y_out <- ego_y_pos$result
neg_y_out <- ego_y_neg$result

# 把所有 list 列（比如 intersection）转成用逗号连接的字符串
pos_y_out[] <- lapply(pos_y_out, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})
neg_y_out[] <- lapply(neg_y_out, function(x) {
  if (is.list(x)) vapply(x, function(z) paste(z, collapse = ","), character(1)) else x
})

write.csv(pos_y_out, "GO_enrichment_yellow_pos_BP.csv", row.names = FALSE)
write.csv(neg_y_out, "GO_enrichment_yellow_neg_BP.csv", row.names = FALSE)
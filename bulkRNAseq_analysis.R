# title: "bulk RNA-seq analysis of Tpam vs. conv. T cells"
# author: "Veronika Cimermanova"
# date: "2023-04-23"

# import packages
library(ggplot2)
library(DESeq2)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(readr)
library(dplyr)
library(circlize)
library(tidyverse)
library(tidyr)
library(cluster)
library(ggsignif)
library(ungeviz)
library(multcomp)
library(ggrepel)
library(RColorBrewer)

#### import of count.table and meta.data ####
# download provided files at:, and create your own path to the files here below
count.table <- read.csv("~/count.table.csv", row.names=1, sep=";")
meta.data <- read.csv("~/meta.data.csv", row.names=1, sep = ";")

# split of data to perform analysis on individual time points
counts_d0 <- dplyr::select(count.table, c(1:6))
counts_d1 <- dplyr::select(count.table, c(7:12))
counts_d3 <- dplyr::select(count.table, c(13:18))

meta.data_d0 <- meta.data[c(1:6), c(1:3)]
meta.data_d1 <- meta.data[c(7:12), c(1:3)]
meta.data_d3 <- meta.data[c(13:18), c(1:3)]

#### DeSeq2 ANALYSIS  on individual time points ####
# DeSeq2 with Tpam genes as upregulated and conv. genes as downregulated
d0_ddsFullCounts <- DESeqDataSetFromMatrix(
  countData = counts_d0,
  colData = meta.data_d0,
  design = ~ Cell_type)
d0_dds <- DESeq(d0_ddsFullCounts)

d1_ddsFullCounts <- DESeqDataSetFromMatrix(
  countData = counts_d1,
  colData = meta.data_d1,
  design = ~ Cell_type)
d1_dds <- DESeq(d1_ddsFullCounts)

d3_ddsFullCounts <- DESeqDataSetFromMatrix(
  countData = counts_d3,
  colData = meta.data_d3,
  design = ~ Cell_type)
d3_dds <- DESeq(d3_ddsFullCounts)

# converting results to data frames, adding FDR(-log10(P-value)) and gene names
d0_R <- as.data.frame(results(d0_dds))
d0_R$FDR <- -log10(d0_R$padj)
d0_R$symbol <- mapIds(org.Mm.eg.db, keys = rownames(d0_R), keytype = "ENSEMBL", column = "SYMBOL")

d1_R <- as.data.frame(results(d1_dds))
d1_R$FDR <- -log10(d1_R$padj)
d1_R$symbol <- mapIds(org.Mm.eg.db, keys = rownames(d1_R), keytype = "ENSEMBL", column = "SYMBOL")

d3_R <- as.data.frame(results(d3_dds))
d3_R$FDR <- -log10(d3_R$padj)
d3_R$symbol <- mapIds(org.Mm.eg.db, keys = rownames(d3_R), keytype = "ENSEMBL", column = "SYMBOL")

## identification of DEGs
# creation of extra column in data frame
d0_R$DEG <- "NO"
d1_R$DEG <- "NO"
d3_R$DEG <- "NO"

# assigning as DEGs only genes which have > 2 (-log10 P-value) and > 1 or < -1 (log2 Fold Change) upregulated/downregulated
d0_R$DEG[d0_R$FDR> 2.0 & d0_R$log2FoldChange> 1.0] <- "UP"
d0_R$DEG[d0_R$FDR> 2.0 & d0_R$log2FoldChange< -1.0] <- "DOWN"

d1_R$DEG[d1_R$FDR> 2.0 & d1_R$log2FoldChange> 1.0] <- "UP"
d1_R$DEG[d1_R$FDR> 2.0 & d1_R$log2FoldChange< -1.0] <- "DOWN"

d3_R$DEG[d3_R$FDR> 2.0 & d3_R$log2FoldChange> 1.0] <- "UP"
d3_R$DEG[d3_R$FDR> 2.0 & d3_R$log2FoldChange< -1.0] <- "DOWN"

# extraction of the TOP 30 DEGs upregulated/downregulated in each time point
# split for upregulated and downregulated
UP_d0 <- subset(d0_R, DEG == "UP")
DOWN_d0 <- subset(d0_R, DEG == "DOWN")

UP_d1 <- subset(d1_R, DEG == "UP")
DOWN_d1 <- subset(d1_R, DEG == "DOWN")

UP_d3 <- subset(d3_R, DEG == "UP")
DOWN_d3 <- subset(d3_R, DEG == "DOWN")

# ordering genes by highest FDR
UP_d0 <- UP_d0[order(UP_d0$FDR, decreasing = TRUE),]
DOWN_d0 <- DOWN_d0[order(DOWN_d0$FDR, decreasing = TRUE),]

UP_d1 <- UP_d1[order(UP_d1$FDR, decreasing = TRUE),]
DOWN_d1 <- DOWN_d1[order(DOWN_d1$FDR, decreasing = TRUE),]

UP_d3 <- UP_d3[order(UP_d3$FDR, decreasing = TRUE),]
DOWN_d3 <- DOWN_d3[order(DOWN_d3$FDR, decreasing = TRUE),]

# extracting first 30 genes into a file
write_csv(UP_d0[c(1:30), c(1:8)], file = "UP_d0.csv")
write_csv(DOWN_d0[c(1:30), c(1:8)], file = "DOWN_d0.csv")

write_csv(UP_d1[c(1:30), c(1:8)], file = "UP_d1.csv")
write_csv(DOWN_d1[c(1:30), c(1:8)], file = "DOWN_d1.csv")

write_csv(UP_d3[c(1:30), c(1:8)], file = "UP_d3.csv")
write_csv(DOWN_d3[c(1:30), c(1:8)], file = "DOWN_d3.csv")

## detecting similar DEGs between groups
# extracting gene names
UP_d0_genes <- UP_d0$symbol
DOWN_d0_genes <- DOWN_d0$symbol
UP_d1_genes <- UP_d1$symbol
DOWN_d1_genes <- DOWN_d1$symbol
UP_d3_genes <- UP_d3$symbol
DOWN_d3_genes <- DOWN_d3$symbol

# comparing each group to other groups
UP_d0_vs_UP_d1 <- intersect(UP_d0_genes, UP_d1_genes)
UP_d0_vs_DOWN_d1 <- intersect(UP_d0_genes, DOWN_d1_genes)
UP_d0_vs_UP_d3 <- intersect(UP_d0_genes, UP_d3_genes)
UP_d0_vs_DOWN_d3 <- intersect(UP_d0_genes, DOWN_d3_genes)
DOWN_d0_vs_UP_d1 <- intersect(DOWN_d0_genes, UP_d1_genes)
DOWN_d0_vs_DOWN_d1 <- intersect(DOWN_d0_genes, DOWN_d1_genes)
DOWN_d0_vs_UP_d3 <- intersect(DOWN_d0_genes, UP_d3_genes)
DOWN_d0_vs_DOWN_d3 <- intersect(DOWN_d0_genes, DOWN_d3_genes)
UP_d1_vs_UP_d3 <- intersect(UP_d1_genes, UP_d3_genes)
UP_d1_vs_DOWN_d3 <- intersect(UP_d1_genes, DOWN_d3_genes)
DOWN_d1_vs_UP_d3 <- intersect(DOWN_d1_genes, UP_d3_genes)
DOWN_d1_vs_DOWN_d3 <- intersect(DOWN_d1_genes, DOWN_d3_genes)
  # results of these comparisons were inspected manually by rewriting it into a table
  # we discovered Ikzf2 as the only upregulated DEG in Tpam cells throughout time points
  # we found that 112 DEGs of Tpam cells on d0 were on d3 increased in conv. cells
  # and that 90 DEGs of conv. cells on d0 were on d3 increased in Tpam cells

#### Volcano plots ####
## d0
# creation of extra column in data frame to choose gene names for volcano plot
d0_R$gene_plot <- "NO"
d0_R$gene_plot[d0_R$FDR> 80.0 & d0_R$log2FoldChange> 1.0] <- "Yes"
d0_R$gene_plot[d0_R$FDR> 60.0 & d0_R$log2FoldChange< -1.0] <- "Yes"
d0_R$gene_plot[d0_R$FDR> 50.0 & d0_R$log2FoldChange< -2.2] <- "Yes"
d0_R$gene_plot[d0_R$FDR> 25.0 & d0_R$log2FoldChange> 5] <- "Yes"

# creation of extra column to rewrite wanted gane names for volcano plot
d0_R$delabel <- NA
d0_R$delabel[d0_R$gene_plot != "NO"] <- d0_R$symbol[d0_R$gene_plot != "NO"]

# Volcano plot
# Thesis = "Figure 10A"
d0_plot <- ggplot(data=d0_R, aes(x=log2FoldChange, y=FDR, col=DEG, label=delabel))+
  geom_point(size=1)+
  geom_label_repel(
    force=5.0,
    alpha=1,
    max.overlaps = Inf,
    fontface="bold",
    size=6,
    min.segment.length = 0,
    box.padding = unit(0.8, "lines"))+
  scale_color_manual(values=c("#2166AC", "#666666", "#B2182B"))+ #colours assigned alphabetically
  geom_vline(xintercept=c(-1.0, 1.0), col="black", lty="dashed")+
  geom_hline(yintercept=2.0, col="black", lty="dashed")+
  theme_light()+
  theme(plot.title = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  xlab("log2(Fold Change)")+
  ylab("-log10(P-value)")+
  theme(legend.position = "none")+
  ggtitle("day 0 (1,988 DEGs)")
print(d0_plot)

## d1
# creation of extra column in data frame to choose gene names for volcano plot
d1_R$gene_plot <- "NO"
d1_R$gene_plot[d1_R$FDR> 7.5 & d1_R$log2FoldChange> 1.2] <- "Yes"
d1_R$gene_plot[d1_R$FDR> 5.0 & d1_R$log2FoldChange< -1.0] <- "Yes"

# creation of extra column to rewrite wanted gane names for volcano plot
d1_R$delabel <- NA
d1_R$delabel[d1_R$gene_plot != "NO"] <- d1_R$symbol[d1_R$gene_plot != "NO"]

# Volcano plot
# Thesis = "Figure 10B"
d1_plot <- ggplot(data=d1_R, aes(x=log2FoldChange, y=FDR, col=DEG, label=delabel))+
  geom_point()+
  geom_label_repel(
    force=5.0,
    alpha=1,
    max.overlaps = Inf,
    fontface="bold",
    size=6,
    min.segment.length = 0,
    box.padding = unit(0.8, "lines"))+
  scale_color_manual(values=c("#2166AC", "#666666", "#B2182B"))+ #colours assigned alphabetically
  geom_vline(xintercept=c(-1.0, 1.0), col="black", lty="dashed")+
  geom_hline(yintercept=2.0, col="black", lty="dashed")+
  theme_light()+
  theme(plot.title = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  xlab("log2(Fold Change)")+
  ylab("-log10(P-value)")+
  theme(legend.position = "none")+
  ggtitle("day 1 (104 DEGs)")
print(d1_plot)

## d3
# creation of extra column in data frame to choose gene names for volcano plot
d3_R$gene_plot <- "NO"
d3_R$gene_plot[d3_R$FDR> 15.0 & d3_R$log2FoldChange> 1.2] <- "Yes"
d3_R$gene_plot[d3_R$FDR> 15.0 & d3_R$log2FoldChange< -1.0] <- "Yes"
d3_R$gene_plot[d3_R$FDR> 5.0 & d3_R$log2FoldChange> 3.5] <- "Yes"
d3_R$gene_plot[d3_R$FDR> 5.0 & d3_R$log2FoldChange< -3.0] <- "Yes"

# creation of extra column to rewrite wanted gane names for volcano plot
d3_R$delabel <- NA
d3_R$delabel[d3_R$gene_plot != "NO"] <- d3_R$symbol[d3_R$gene_plot != "NO"]

# Volcano plot
# Thesis = "Figure 10B"
d3_plot <- ggplot(data=d3_R, aes(x=log2FoldChange, y=FDR, col=DEG, label=delabel))+
  geom_point()+
  geom_label_repel(
    force=5.0,
    alpha=1,
    max.overlaps = Inf,
    fontface="bold",
    size=6,
    min.segment.length = 0,
    box.padding = unit(0.8, "lines"))+
  scale_color_manual(values=c("#2166AC", "#666666", "#B2182B"))+ #colours assigned alphabetically
  geom_vline(xintercept=c(-1.0, 1.0), col="black", lty="dashed")+
  geom_hline(yintercept=2.0, col="black", lty="dashed")+
  theme_light()+
  theme(plot.title = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  xlab("log2(Fold Change)")+
  ylab("-log10(P-value)")+
  theme(legend.position = "none")+
  ggtitle("day 3 (309 DEGs)")
print(d3_plot)

#### Heatmap of the steady-state signature Tpam genes ####
# scRNA-seq DEGs of Tpam cluster
genes_Tpam_cluster <- c("Nfkbid", "Tnf", "Myc", "Dusp2", "Nab2", "Tnfrsf9", "Txnip", "Klf2",
                        "Cd69", "Ier2","Tsc22d3", "Eif4a1","Klf3", "Zfp36", "Ccl5", "Srm", "Il7r",
                        "Kdm6b", "Rilpl2","S100a10", "Zfp36l1", "Ptpn6", "Relb", "Orai1", "Nfkbia",
                        "Srsf2", "Mat2a", "Odc1","Gnl3", "Hspa5", "S100a6", "Ifngr1", "Bcl3", "Ppan",
                        "Ctsd", "Nolc1","Ranbp1", "Cxcr3", "Nhp2", "Cd7", "Pa2g4", "Sh2d1a", "Ahnak",
                        "Hmgb2", "Tagap","Klf6", "Ifitm10", "Prr13")
# extracting rownames of data frame based on scRNA-seq genes of Tpam cluster
d0_genes <- d0_R[d0_R$symbol %in% genes_Tpam_cluster, ]
d1_genes <- d1_R[d1_R$symbol %in% genes_Tpam_cluster, ]
d3_genes <- d3_R[d3_R$symbol %in% genes_Tpam_cluster, ]

# making dds counts for scRNA-seq genes of Tpam cluster
d0_dds_matrix <- d0_dds[row.names(d0_dds) %in% row.names(d0_genes)]
d1_dds_matrix <- d1_dds[row.names(d1_dds) %in% row.names(d1_genes)]
d3_dds_matrix <- d3_dds[row.names(d3_dds) %in% row.names(d3_genes)]

# normalization of counts for Z value and adding names to columns based on meta.data
d0_map_matrix <- counts(d0_dds_matrix, normalized = T)
d0_map_matrix_Z <- t(apply(d0_map_matrix, 1, scale))
colnames(d0_map_matrix_Z) <- rownames(meta.data_d0)

d1_map_matrix <- counts(d1_dds_matrix, normalized = T)
d1_map_matrix_Z <- t(apply(d1_map_matrix, 1, scale))
colnames(d1_map_matrix_Z) <- rownames(meta.data_d1)

d3_map_matrix <- counts(d3_dds_matrix, normalized = T)
d3_map_matrix_Z <- t(apply(d3_map_matrix, 1, scale))
colnames(d3_map_matrix_Z) <- rownames(meta.data_d3)

# Heatmaps
# Thesis = "Figure 10D"
col_names_timepoints <- c(" ",  " ", "Tpam"," ", " ", "conv.")

H0 <- Heatmap(d0_map_matrix_Z, cluster_rows = T, 
              cluster_columns = F, col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
              column_labels = col_names_timepoints,
              column_names_rot = 0,
              column_names_gp = gpar(fontsize = 20),
              row_labels = d0_R[rownames(d0_genes),]$symbol,
              row_names_gp = gpar(fontsize = 16, fontface = "italic"),
              column_title = "day 0",
              column_title_gp = gpar(fontsize = 20),
              column_dend_side = "bottom", row_dend_width = unit(2, "cm"),
              name = "Z score")

H1 <- Heatmap(d1_map_matrix_Z, cluster_rows = T, 
              cluster_columns = F, col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
              column_labels = col_names_timepoints,
              column_names_rot = 0,
              column_names_gp = gpar(fontsize = 20),
              row_labels = d1_R[rownames(d1_genes),]$symbol,
              row_names_gp = gpar(fontsize = 16, fontface = "italic"),
              column_title = "day 1",
              column_title_gp = gpar(fontsize = 20),
              column_dend_side = "bottom",  row_dend_width = unit(2, "cm"),
              name = "Z score")

H3 <- Heatmap(d3_map_matrix_Z, cluster_rows = T, 
              cluster_columns = F, col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
              column_labels = col_names_timepoints,
              column_names_rot = 0,
              column_names_gp = gpar(fontsize = 20),
              row_labels = d3_R[rownames(d3_genes),]$symbol,
              row_names_gp = gpar(fontsize = 16, fontface = "italic"),
              column_title = "day 3",
              column_title_gp = gpar(fontsize = 20),
              column_dend_side = "bottom", row_dend_width = unit(2, "cm"),
              name = "Z score")

# final combined heatmaps
H0 + H1 + H3

#### RNA expression in time for Mki67 and Ikzf2 ####
## DeSeq2 analysis throughout time
ddsFullCounts <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = meta.data,
  design = ~ Cell_type)
dds <- DESeq(ddsFullCounts)

# converting results to data frames, adding FDR(-log10(P-value)) and gene names
R <- as.data.frame(results(dds))
R$FDR <- -log10(R$padj)
R$symbol <- mapIds(org.Mm.eg.db, keys = rownames(R), keytype = "ENSEMBL", column = "SYMBOL")

# defying wanted genes
Mki67  <- subset(R, symbol == "Mki67")
Ikzf2 <- subset(R, symbol == "Ikzf2")

# making dds counts similar length for genes
Mki67_dds_matrix <- dds[row.names(dds) %in% row.names(Mki67)]
Ikzf2_dds_matrix <- dds[row.names(dds) %in% row.names(Ikzf2)]

# log2 normalization
Mki67_log2_NM <- t(log2((counts(Mki67_dds_matrix[rownames(Mki67), ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(Mki67_dds_matrix), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(rownames(Mki67))+1):ncol(.))
Mki67_log2_NM$symbol <- mapIds(org.Mm.eg.db, keys = Mki67_log2_NM$gene, keytype = "ENSEMBL", column = "SYMBOL")

Ikzf2_log2_NM <- t(log2((counts(Ikzf2_dds_matrix[rownames(Ikzf2), ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(Ikzf2_dds_matrix), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(rownames(Ikzf2))+1):ncol(.))
Ikzf2_log2_NM$symbol <- mapIds(org.Mm.eg.db, keys = Ikzf2_log2_NM$gene, keytype = "ENSEMBL", column = "SYMBOL")

# changing order for graphs to have Tpam cells first
Mki67_log2_NM$Cell_type <- factor(Mki67_log2_NM$Cell_type, levels = rev(levels(Mki67_log2_NM$Cell_type)))
Ikzf2_log2_NM$Cell_type <- factor(Ikzf2_log2_NM$Cell_type, levels = rev(levels(Ikzf2_log2_NM$Cell_type)))

### GRAPHS creation and extra statistic ###
## Mki67
# Thesis = "Figure 8D"
# making summary statistics
sum_Mki67 <-group_by(Mki67_log2_NM, gene, symbol, Time, Cell_type) %>%
  summarise(
    count = n(),
    mean = mean(expression, na.rm = TRUE),
    sd = sd(expression, na.rm = TRUE))
sum_Mki67$expression <- sum_Mki67$mean

# statistic for separate time expression for graphs
d0_Mki67 <- subset(Mki67_log2_NM, Time == "0")
d1_Mki67 <- subset(Mki67_log2_NM, Time == "1")
d3_Mki67 <- subset(Mki67_log2_NM, Time == "3")

# check of distribution
with(d0_Mki67, shapiro.test(expression[Cell_type == "Tpam"]))
with(d0_Mki67, shapiro.test(expression[Cell_type == "Conv"]))
with(d1_Mki67, shapiro.test(expression[Cell_type == "Tpam"]))
with(d1_Mki67, shapiro.test(expression[Cell_type == "Conv"]))
with(d3_Mki67, shapiro.test(expression[Cell_type == "Tpam"]))
with(d3_Mki67, shapiro.test(expression[Cell_type == "Conv"]))
# result: normal distribution

#Unpaired T-test
t.test(expression ~ Cell_type, data = d0_Mki67, var.equal = TRUE)
t.test(expression ~ Cell_type, data = d1_Mki67, var.equal = TRUE) # 0.0048
t.test(expression ~ Cell_type, data = d3_Mki67, var.equal = TRUE)

# plot
Mki67_plot <- ggplot(data=Mki67_log2_NM, aes(x=Time, y=expression, color=Cell_type))+
    geom_point(aes(color=Cell_type), position=position_dodge(width=1), size=4, alpha=0.7)+
    geom_line(sum_Mki67, mapping=aes(x=Time, y=expression, color=Cell_type),
              position=position_dodge(width=1), size=0.75, linetype = "dashed")+
    theme_classic()+
    ylim(0,15)+
    xlab("Time point (day)")+
    ylab("Mki67 RNA expression\n(log2 normalized counts)")+
    labs(shape = "Cells")+
    scale_color_manual(values=c("forestgreen", "gray20"), labels=c("Tpam" = "Tpam", "Conv" = "conv."))+
    theme(axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.x = element_text(size=20, color="black"),
          axis.text.y = element_text(size=20, color="black"),
          legend.title = element_text(size=20,),
          legend.text = element_text(size=20))
Mki67_plot <- Mki67_plot + geom_hpline(sum_Mki67, mapping=aes(x=Time, y=expression, color=Cell_type),
                                       width=0.3, position=position_dodge(width=1))
Mki67_plot <- Mki67_plot + geom_signif(y_position = c(11,13.5,14), xmin = c(-0.25,0.75, 2.75), 
                           xmax = c(0.25,1.25, 3.25), annotation = c("NS", "p=0.048", "NS"),
                           tip_length = 0.03, textsize = 5)
print(Mki67_plot)

## Ikzf2
# Thesis = "Figure 10C"
# making summary statistics
sum_Ikzf2 <-group_by(Ikzf2_log2_NM, gene, symbol, Time, Cell_type) %>%
  summarise(
    count = n(),
    mean = mean(expression, na.rm = TRUE),
    sd = sd(expression, na.rm = TRUE))
sum_Ikzf2$expression <- sum_Ikzf2$mean

# statistic for separate time expression for graphs
d0_Ikzf2 <- subset(Ikzf2_log2_NM, Time == "0")
d1_Ikzf2 <- subset(Ikzf2_log2_NM, Time == "1")
d3_Ikzf2 <- subset(Ikzf2_log2_NM, Time == "3")

# check of distribution
with(d0_Ikzf2, shapiro.test(expression[Cell_type == "Tpam"]))
with(d0_Ikzf2, shapiro.test(expression[Cell_type == "Conv"]))
with(d1_Ikzf2, shapiro.test(expression[Cell_type == "Tpam"]))
with(d1_Ikzf2, shapiro.test(expression[Cell_type == "Conv"]))
with(d3_Ikzf2, shapiro.test(expression[Cell_type == "Tpam"]))
with(d3_Ikzf2, shapiro.test(expression[Cell_type == "Conv"]))
# result: normal distribution

#Unpaired T-test
t.test(expression ~ Cell_type, data = d0_Ikzf2, var.equal = TRUE) # 0.005
t.test(expression ~ Cell_type, data = d1_Ikzf2, var.equal = TRUE) # 0.0084
t.test(expression ~ Cell_type, data = d3_Ikzf2, var.equal = TRUE) # 0.0011

# plot
Ikzf2_plot <- ggplot(data=Ikzf2_log2_NM, aes(x=Time, y=expression, color=Cell_type))+
    geom_point(aes(color=Cell_type), position=position_dodge(width=1), size=4, alpha=0.7)+
    geom_line(sum_Ikzf2, mapping=aes(x=Time, y=expression, color=Cell_type),
              position=position_dodge(width=1), size=0.75, linetype = "dashed")+
    theme_classic()+
    ylim(0,15)+
    xlab("Time point (day)")+
    ylab("Ikzf2 RNA expression\n(log2 normalized counts)")+
    labs(shape = "Cells")+
    scale_color_manual(values=c("forestgreen", "gray20"), labels=c("Tpam" = "Tpam", "Conv" = "conv."))+
    theme(axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.x = element_text(size=20, color="black"),
          axis.text.y = element_text(size=20, color="black"),
          legend.title = element_text(size=20,),
          legend.text = element_text(size=20))
Ikzf2_plot <- Ikzf2_plot + geom_hpline(sum_Ikzf2, mapping=aes(x=Time, y=expression, color=Cell_type),
                                       width=0.3, position=position_dodge(width=1))
Ikzf2_plot <- Ikzf2_plot + geom_signif(y_position = c(10,12,14), xmin = c(-0.25,0.75, 2.75), 
                           xmax = c(0.25,1.25, 3.25), annotation = c("p=0.005", "p=0.0084", "p=0.0011"),
                           tip_length = 0.03, textsize = 5)
print(Ikzf2_plot)

#### session info ####
sessionInfo()
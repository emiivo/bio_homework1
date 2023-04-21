#!/bin/bash

RESULTS_DIR=~/HW1/results
OUTPUT_DIR=~/HW1/outputs

Rscript -e '

library(DESeq2)

library(ggplot2)

Counts <- read.delim("~/HW1/results/correct_non_mapping_counts.tsv", header = TRUE, row.names = NULL, sep = "\t")

Counts <- Counts[,-1]

colnames(Counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")

Counts <- round(Counts) # round Counts to nearest integer

Counts <- Counts[which(rowSums(Counts) > 50),]

condition <- factor(c("1", "2", "3", "4"))

coldata <- data.frame(row.names = colnames(Counts), condition)
#Counts <- as.numeric(unlist(Counts))
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

dds <- DESeq(dds)

vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata, intgroup = "condition")

DEGs <- results(dds)

if ("gene_symbols" %in% colnames(DEGs)) {
  DEGs <- DEGs[, c("gene_symbols", "log2FoldChange", "lfcSE", "pvalue", "padj")]
} else {
  DEGs <- DEGs[, c("log2FoldChange", "lfcSE", "pvalue", "padj")]
}

DEGs <- DEGs[which(DEGs$padj < 0.05),]

DEGs <- DEGs[order(abs(DEGs$log2FoldChange), decreasing = TRUE),]

write.table(DEGs, file = "mapping_based_DEGs.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

DEGs <- read.delim("mapping_based_DEGs.tsv", header = TRUE, row.names = 1, sep = "\t")

# Define cutoff values for significance and log2 fold change
pval_cutoff <- 0.05
log2fc_cutoff <- 1

# Create the volcano plot

volcano_plot <- ggplot(data = DEGs, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= log2fc_cutoff & padj <= pval_cutoff, "DEG", "Not DEG")), size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("DEG" = "red", "Not DEG" = "black")) +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "gray", size = 0.5) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "gray", size = 0.5) +
  scale_x_continuous(name = "Log2 Fold Change", limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
  scale_y_continuous(name = "-log10(Adjusted P-value)", limits = c(0, 50), breaks = seq(0, 50, 10)) +
  theme_bw()

# Save the plot as a PNG file

ggsave("volcano_plot.png", plot = volcano_plot, width = 7, height = 5, dpi = 300)

DEGs <- read.table("'${RESULTS_DIR}/mapping_based_DEGs.tsv'", header = TRUE, sep = "\t", row.names = 1)

# Select top 50 DEGs by adjusted p-value
top_DEGs <- DEGs[order(DEGs$padj)[1:50], ]

# Create heatmap of selected DEGs
heatmap.2(as.matrix(top_DEGs), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          scale = "row",
          margins = c(12, 9),
          labCol = colnames(top_DEGs),
          labRow = rownames(top_DEGs),
          col = colorRampPalette(c("navyblue", "white", "firebrick2"))(100))'
          

Rscript -e '

library(DESeq2)

library(ggplot2)


Counts <- read.delim("~/HW1/outputs/correct_gene_count_matrix.csv", header = TRUE, row.names = 1, sep = ",")

Counts <- Counts[which(rowSums(Counts) > 50),]

condition <- factor(c("1", "2", "3", "4"))

coldata <- data.frame(row.names = colnames(Counts), condition)

dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

dds <- DESeq(dds)

vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata, intgroup = "condition")

DEGs <- results(dds)

if ("gene_symbols" %in% colnames(DEGs)) {
  DEGs <- DEGs[, c("gene_symbols", "log2FoldChange", "lfcSE", "pvalue", "padj")]
} else {
  DEGs <- DEGs[, c("log2FoldChange", "lfcSE", "pvalue", "padj")]
}

DEGs <- DEGs[which(DEGs$padj < 0.05),]

DEGs <- DEGs[order(abs(DEGs$log2FoldChange), decreasing = TRUE),]

write.table(DEGs, file = "alignment_free_DEGs.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

DEGs <- read.delim("alignment_free_DEGs.tsv", header = TRUE, row.names = 1, sep = "\t")

# Define cutoff values for significance and log2 fold change
pval_cutoff <- 0.05
log2fc_cutoff <- 1

# Create the volcano plot

volcano_plot <- ggplot(data = DEGs, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= log2fc_cutoff & padj <= pval_cutoff, "DEG", "Not DEG")), size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("DEG" = "red", "Not DEG" = "black")) +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "gray", size = 0.5) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "gray", size = 0.5) +
  scale_x_continuous(name = "Log2 Fold Change", limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
  scale_y_continuous(name = "-log10(Adjusted P-value)", limits = c(0, 50), breaks = seq(0, 50, 10)) +
  theme_bw()

# Save the plot as a PNG file

ggsave("volcano_plot.png", plot = volcano_plot, width = 7, height = 5, dpi = 300)

DEGs <- read.table("'${RESULTS_DIR}/alignment_free_DEGs.tsv'", header = TRUE, sep = "\t", row.names = 1)

# Select top 50 DEGs by adjusted p-value
top_DEGs <- DEGs[order(DEGs$padj)[1:50], ]

# Create heatmap
heatmap.2(as.matrix(top_DEGs), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          scale = "row",
          margins = c(12, 9),
          labCol = colnames(top_DEGs),
          labRow = rownames(top_DEGs),
          col = colorRampPalette(c("navyblue", "white", "firebrick2"))(100))'
          



Rscript -e '

library(DESeq2)
library(ggplot2)

AF_DEGs <- read.table("'${OUTPUT_DIR}/alignment_free_DEGs.tsv'", header = TRUE, sep = "\t", row.names = 1)
MB_DEGs <- read.table("'${RESULTS_DIR}/mapping_based_DEGs.tsv'", header = TRUE, sep = "\t", row.names = 1)

# Get common DEGs
common_DEGs <- intersect(row.names(AF_DEGs), row.names(MB_DEGs))

# Get DEGs specific to each method
AF_specific_DEGs <- setdiff(row.names(AF_DEGs), common_DEGs)
MB_specific_DEGs <- setdiff(row.names(MB_DEGs), common_DEGs)

# Create Venn diagram
venn.diagram(list(AF = row.names(AF_DEGs), MB = row.names(MB_DEGs)),
             filename = "'${RESULTS_DIR}/DEG_venn_diagram.png'",
             category.names = c("Alignment-free", "Mapping-based"),
             fill = c("cornflowerblue", "darkorange1"),
             alpha = c(0.8, 0.8),
             lty = "blank",
             cex = 1.5,
             fontface = "bold",
             cat.cex = 1.5,
             cat.fontface = "bold",
             euler.d = TRUE,
             count = TRUE,
             circle.col = c("cornflowerblue", "darkorange1"),
             print.mode = "single",
             main = "DEGs Venn Diagram")

# Add labels for common and specific DEGs
grid.text(label = paste("Common DEGs:\n", length(common_DEGs)), x = 0.58, y = 0.5, gp = gpar(fontsize = 12, fontface = "bold", col = "black"))
grid.text(label = paste("Alignment-free specific DEGs:\n", length(AF_specific_DEGs)), x = 0.25, y = 0.25, gp = gpar(fontsize = 12, fontface = "bold", col = "cornflowerblue"))
grid.text(label = paste("Mapping-based specific DEGs:\n", length(MB_specific_DEGs)), x = 0.91, y = 0.25, gp = gpar(fontsize = 12, fontface = "bold", col = "darkorange1"))'

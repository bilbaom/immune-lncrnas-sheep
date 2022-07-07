library(RColorBrewer)
library("viridis")
setwd("/home/labo/datuak/atlas/")
library("DESeq2")
library("ggplot2")
library(factoextra)
library("Rtsne")
library("apeglm")
library("gprofiler2")

# metadata table
colData <- read.csv("immunetissues_exp.csv")
colData <- colData[!apply(is.na(colData) | colData == "", 1, all),]
colData <- colData[c("Experiment", "Status", "Tissue", "Tissue2", "Age", "Breed", "Sex", "Library", "BioProject", "Condition", "Effect")]
rownames(colData) <- colData$Experiment
#write.table(colData,"coldata.csv", sep = ",")

# Open saved expression matrix of raw counts (exported from swichanalyser)
counts <- read.csv("expr/rawcounts.csv", row.names = 1)

# keep only ENS genes and novel lncRNAs
lncids <- read.table("custom_lnc/lnc_ids_genelevel_filt.txt", sep = "\t", header=TRUE)
enscounts <- counts[startsWith(rownames(counts), "ENS"),]
lnccounts <- merge(counts, lncids, by.x=0, by.y="ID")
rownames(lnccounts) <- lnccounts$Row.names
lnccounts <- lnccounts[,-1]
counts = rbind(enscounts,lnccounts)
rownames(counts) <- gsub( " .*", "", rownames(counts))
rownames(counts)

# remove rRNA, tRNA, snRNA, snoRNA, mtRNA
rnas <- read.csv("references/rrna_trna_snrna.txt", row.names = 1, sep="\t")
counts <- counts[!(rownames(counts) %in% rownames(rnas)),]

### ALL SAMPLES #######################################################################

# Filter reads
counts_filtered <- counts[rowMeans(counts)>10 ,]

# Differential expression DESeq2

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered), 
  colData = colData, 
  design = ~ Tissue2 + Age + Sex + Library + Status)

dds <- DESeq(dds)
res <- results(dds)
normalized_counts <- counts(dds, normalized=TRUE)

contrast <- c("Status", "Treatment", "Control")

# unshrunken FC
res_unshrunken <- results(dds, contrast=contrast, alpha = 0.05)
summary(res_unshrunken)
r = as.data.frame(res_unshrunken)
r = r[r$padj < 0.05,] # adjusted p value (FDR)
r = r[abs(r$log2FoldChange) > 0.5849 ,] # logFC
r = na.omit(r)

write.table(r,"expr/de_results_p05_fc58_.txt", sep = "\t")

# shrunken FC
res_shrunken <- lfcShrink(dds, coef="Status_Treatment_vs_Control", res=res_unshrunken)
summary(res_shrunken)
r2 = as.data.frame(res_shrunken)
r2 = r2[r2$padj < 0.05,] # adjusted p value (FDR)
r2 = r2[abs(r2$log2FoldChange) > 0.5849 ,] # logFC
r2 = na.omit(r2)

### BLOOD SAMPLES #######################################################################
# Filter metadata 
colData_b <- colData[colData["Tissue"]=="Blood",]
colData_b <- rbind(colData_b, colData[colData["Tissue"]=="Cells",]) # ADD CELL SAMPLES
colData_b$Tiss_BP <- paste(colData_b$Tissue2,colData_b$BioProject)

# Remove sample groups without effect in publication
colData_b <- colData_b[colData_b["Effect"]!="No",]

# Filter reads
counts_b <- counts[,rownames(colData_b) ]
#counts_filtered_b <- counts_b[rowMeans(counts_b)>10 ,]
#counts_filtered_b <- counts_filtered_b[rowSums(counts_filtered_b != 0) > length(counts_b) / 2, ]
counts_filtered_b <- counts_b[rowSums(counts_b > 1) > (length(counts_b) / 2), ]

rownames(colData_b)
colnames(counts_b)

# Differential expression DESeq2

dds_b <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered_b), 
  colData = colData_b, 
  design = ~ Tiss_BP + Status)

dds_b <- estimateSizeFactors(dds_b)
dds_b <- DESeq(dds_b)
nrow(dds_b)
dds_b <- dds_b[which(mcols(dds_b)$betaConv),]
normalized_counts_b <- counts(dds_b, normalized=TRUE)
#write.csv(normalized_counts_b, "expr/normcounts_blood.csv", quote=FALSE)

contrast <- c("Status", "Treatment", "Control")
LOG2FC = log2(1.25)
# unshrunken FC
res_unshrunken_b <- results(dds_b, contrast=contrast, alpha = 0.05)
summary(res_unshrunken_b)
r = as.data.frame(res_unshrunken_b)
r = r[r$padj < 0.05,] # adjusted p value (FDR)
r = r[abs(r$log2FoldChange) > LOG2FC ,] # logFC
r = na.omit(r)

write.table(r,"expr/de_results_p05_fc32_blood.txt", sep = "\t", quote=FALSE)

# shrunken FC
res_shrunken_b <- lfcShrink(dds_b, coef="Status_Treatment_vs_Control", res=res_unshrunken_b)
summary(res_shrunken_b)
r2 = as.data.frame(res_shrunken_b)
r2 = r2[r2$padj < 0.05,] # adjusted p value (FDR)
r2 = r2[abs(r2$log2FoldChange) > LOG2FC ,] # logFC
r2 = na.omit(r2)

write.table(r2,"expr/de_results_p05_fc32_blood_shrunken.txt", sep = "\t", quote=FALSE)

# GO enrichment
r2 <- r2[order(r2$pvalue), ]
r2 <- r2[order(r2$padj), ]
gostres <- gost(query = rownames(r2), organism = "oarambouillet", correction_method = "fdr",
                domain_scope = "custom_annotated", custom_bg = rownames(normalized_counts_b), as_short_link = FALSE)
gores <- apply(gostres$result,2,as.character)
gores <- gores[, colnames(gores) != "parents"]
write.table(gores, "expr/gProfiler_blood.csv", sep = "\t", quote = FALSE, row.names = FALSE)
gostplot(gostres, capped = FALSE, interactive = TRUE)

gnames <- gconvert(query = rownames(r), organism = "oarambouillet",
                   target="ENSG", mthreshold = Inf, filter_na = FALSE)
?gost
### TISSUE SAMPLES #######################################################################
# Filter metadata ((KEEP ONLY LYMPH NODE SAMPLES, THERE IS NOT ANY STIMULATED SPLEEN, THYMUS...))
colData_t <- colData[colData["Tissue"]=="Tissue",]
colData_t <- colData[grep("lymph", colData$Tissue2), ] # ONLY LYMPH NODES
colData_t$Tiss_BP <- paste(colData_t$Tissue2,colData_t$BioProject)

# Filter reads
counts_t <- counts[,rownames(colData_t) ]
#counts_filtered_t <- counts_t[rowMeans(counts_t)>10 ,]
counts_filtered_t <- counts_t[rowSums(counts_t > 1) > (length(counts_t) / 2), ]

# Differential expression DESeq2
dds_t <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered_t), 
  colData = colData_t, 
  design = ~ Tissue2 + Age + Sex + Library + Status)

dds_t <- DESeq(dds_t)
nrow(dds_t)
dds_t <- dds_t[which(mcols(dds_t)$betaConv),]
nrow(dds_t)
normalized_counts_t <- counts(dds_t, normalized=TRUE)
#write.csv(normalized_counts_t, "expr/normcounts_lymph.csv", quote=FALSE)

contrast <- c("Status", "Treatment", "Control")

# unshrunken FC
res_unshrunken_t <- results(dds_t, contrast=contrast, alpha = 0.05)
summary(res_unshrunken_t)
r = as.data.frame(res_unshrunken_t)
r = r[r$padj < 0.05,] # adjusted p value (FDR)
r = r[abs(r$log2FoldChange) > LOG2FC ,] # logFC
r = na.omit(r)

write.table(r,"expr/de_results_p05_fc32_lymph.txt", sep = "\t", quote=FALSE)


# shrunken FC
res_shrunken_t <- lfcShrink(dds_t, coef="Status_Treatment_vs_Control", res=res_unshrunken_t)
summary(res_shrunken_t)
r2 = as.data.frame(res_shrunken_t)
r2 = r2[r2$padj < 0.05,] # adjusted p value (FDR)
r2 = r2[abs(r2$log2FoldChange) > LOG2FC ,] # logFC
r2 = na.omit(r2)

write.table(r2,"expr/de_results_p05_fc32_lymph_shrunken.txt", sep = "\t")

# GO enrichment
r2 <- r2[order(r2$pvalue), ]
r2 <- r2[order(r2$padj), ]
gostres <- gost(query = rownames(r2), organism = "oarambouillet", correction_method = "fdr",
                domain_scope = "custom", custom_bg = rownames(normalized_counts_t), as_short_link = FALSE)
gores <- apply(gostres$result,2,as.character)
gores <- gores[, colnames(gores) != "parents"]
write.table(gores, "expr/gProfiler_lymph.csv", sep = "\t", quote = FALSE, row.names = FALSE)
gostplot(gostres, capped = FALSE, interactive = TRUE)

gnames <- gconvert(query = rownames(r2), organism = "oarambouillet",
                   target="ENSG", mthreshold = Inf, filter_na = FALSE)

### Volcano plot and MA plot ################################################################
### Blood
# Add raw DE results (shrunken or not) and add a column with gene type
de <- as.data.frame(res_shrunken_b)
de$Type <- "No DE"
de$Type[startsWith(rownames(de), "M") & abs(de$log2FoldChange) > LOG2FC & de$padj < 0.05] <- "lncRNA"
de$Type[startsWith(rownames(de), "E") & abs(de$log2FoldChange) > LOG2FC & de$padj < 0.05] <- "Ensembl"

p <- ggplot(data = de, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=Type)) +
  geom_point(alpha=0.4, size=3) +
  scale_color_manual(values=c("tomato", "royalblue","grey"))+
  xlim(c(-2, 2)) +
  geom_vline(xintercept=c(-LOG2FC,LOG2FC),lty=4,col="black",lwd=0.7) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.7) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Blood")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf(file = "expr/volcano_blood.pdf", width = 7, height = 5)
p
dev.off()

# MA plot
plotMA(res_unshrunken_b)
plotMA(res_shrunken_b)

# Add raw DE results (shrunken or not) and add a column with gene type
de <- as.data.frame(res_shrunken_b)
de$Type <- "No DE"
de$Type[startsWith(rownames(de), "M") & abs(de$log2FoldChange) > 0 & de$padj < 0.05] <- "lncRNA"
de$Type[startsWith(rownames(de), "E") & abs(de$log2FoldChange) > 0 & de$padj < 0.05] <- "Ensembl"

ma <- ggplot(data = de, 
             aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 colour=Type)) +
  geom_point(alpha=0.8, size=1) +
  scale_color_manual(values=c("tomato", "royalblue","grey"))+
  geom_hline(yintercept = LOG2FC ,lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -LOG2FC ,lty=4,col="black",lwd=0.5) +
  labs(x="log10 (mean expression)",
       y="log2 (fold change)",
       title="Blood")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf(file = "expr/ma_blood.pdf", width = 7, height = 5)
ma
dev.off()

### Lymph nodes
# Add raw DE results (shrunken or not)
de <- as.data.frame(res_shrunken_t)

# Add a column with gene type
de$Type <- "No DE"
de$Type[startsWith(rownames(de), "M") & abs(de$log2FoldChange) > LOG2FC & de$padj < 0.05] <- "lncRNA"
de$Type[startsWith(rownames(de), "E") & abs(de$log2FoldChange) > LOG2FC & de$padj < 0.05] <- "Ensembl"

p <- ggplot(data = de, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=Type)) +
  geom_point(alpha=0.4, size=3) +
  scale_color_manual(values=c("tomato", "royalblue","grey"))+
  xlim(c(-4.5, 4.5)) +
  ylim(c(0, 6)) +
  geom_vline(xintercept=c(-LOG2FC,LOG2FC),lty=4,col="black",lwd=0.7) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.7) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Lymph nodes")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf(file = "expr/volcano_lymph.pdf", width = 7, height = 5)
p
dev.off()

# MA plot
plotMA(res_unshrunken_t)
plotMA(res_shrunken_t)

# Add raw DE results (shrunken or not) and add a column with gene type
de <- as.data.frame(res_shrunken_t)
de$Type <- "No DE"
de$Type[startsWith(rownames(de), "M") & abs(de$log2FoldChange) > 0 & de$padj < 0.05] <- "lncRNA"
de$Type[startsWith(rownames(de), "E") & abs(de$log2FoldChange) > 0 & de$padj < 0.05] <- "Ensembl"

ma <- ggplot(data = de, 
             aes(x = log10(baseMean), 
                 y = log2FoldChange, 
                 colour=Type)) +
  geom_point(alpha=0.8, size=1) +
  scale_color_manual(values=c("tomato", "royalblue","grey"))+
  geom_hline(yintercept = LOG2FC ,lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -LOG2FC ,lty=4,col="black",lwd=0.5) +
  labs(x="log10 (mean expression)",
       y="log2 (fold change)",
       title="Lymph nodes")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf(file = "expr/ma_lymph.pdf", width = 7, height = 5)
ma
dev.off()


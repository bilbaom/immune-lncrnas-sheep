library("ggplot2")
library(factoextra)
library(RColorBrewer)
library("viridis")
library("limma")
library("sva")
library(GWENA)
library(magrittr)
library(dcanr)
library(igraph)
library("edgeR")
library("minet")
library(pheatmap)
setwd("/home/labo/datuak/atlas/")

# load metadata tables
colData <- read.csv("immunetissues_exp.csv")
colData <- colData[!apply(is.na(colData) | colData == "", 1, all),]
colData <- colData[c("Experiment", "Status", "Tissue", "Tissue2", "Age", "Breed", "Sex", "Library",
                     "BioProject", "Condition", "Effect")]
rownames(colData) <- colData$Experiment
colData$TissueBP <- paste(colData$Tissue2, colData$BioProject, sep = "_")
colData$TissueBPstatus <- paste(colData$TissueBP, colData$Status, sep = "_")

colData_b <- colData[colData["Tissue"]=="Blood",]
colData_b <- rbind(colData_b, colData[colData["Tissue"]=="Cells",]) # ADD CELL SAMPLES
colData_b <- colData_b[colData_b["Effect"]!="No",] # Remove treatments with no effect in article
colData_b <- colData_b[colData_b["Condition"]!="adjuvant",] #outlier
colData_b <- colData_b[colData_b["Experiment"]!="SRX4017443",] #outlier
#write.table(colData_b,"coex/blood/coldata_blood_coex.csv", sep = ",")
table(colData_b$TissueBPstatus)
table(colData_b$Condition)
table(colData_b$Status)

colData_t <- colData[colData["Tissue"]=="Tissue",]
colData_t <- colData[grep("lymph", colData$Tissue2), ] # ONLY LYMPH NODES
#write.table(colData_t,"coex/lymph/coldata_lymph_coex.csv", sep = ",")
table(colData_t$TissueBPstatus)
table(colData_t$Condition)
table(colData_t$Status)

# load Deseq2 normalised count matrixes
counts_b <- read.csv("expr/normcounts_blood.csv", row.names = 1)
counts_t <- read.csv("expr/normcounts_lymph.csv", row.names = 1)

# filter low expression
counts_b <- counts_b[rowSums(counts_b>1)>=(length(counts_b)/2) ,]
counts_t <- counts_t[rowSums(counts_t>1)>=(length(counts_t)/2) ,]

# Open saved expression matrix of raw counts (exported from swichanalyser)
counts <- read.csv("expr/rawcounts.csv", row.names = 1)

# keep only ENS genes and novel lncRNAs
lncids <- read.table("custom_lnc/lnc_ids_genelevel.txt", sep = "\t", header=TRUE)
enscounts <- counts[startsWith(rownames(counts), "ENS"),]
lnccounts <- merge(counts, lncids, by.x=0, by.y="ID")
rownames(lnccounts) <- lnccounts$Row.names
lnccounts <- lnccounts[,-1]
counts = rbind(enscounts,lnccounts)
rownames(counts) <- gsub( " .*", "", rownames(counts))

# remove rRNA, tRNA, snRNA, snoRNA, mtRNA
rnas <- read.csv("references/rrna_trna_snrna.txt", row.names = 1, sep="\t")
counts <- counts[!(rownames(counts) %in% rownames(rnas)),]

########################################################################################
### BLOOD ##############################################################################
########################################################################################
load("coex/coex_atlas.RData")
## Test two normalisation methods: 1: Deseq2 VST 2: EdgeR CTF + asinh

# Design for covariate correction
design_b <- model.matrix(~Status + TissueBP, data=colData_b)
treatment.design_b <- design_b[,1:2]
batch.design_b <- design_b[,-(1:2)]

## EdgeR CTF + asinh
# CTF + asinh normalisation
rawcounts_b <- counts[,rownames(colData_b) ]
rawcounts_b <- rawcounts_b[rownames(counts_b), ]
norm_factors_b <- calcNormFactors(object = rawcounts_b, lib.size = colSums(rawcounts_b), method = "TMM")
counts_b_CTF <- sweep(rawcounts_b, 2, norm_factors_b, "/")
counts_b_CTF <- asinh(counts_b_CTF)
# Correct counts for variation not interesting
counts_b_CTF_corr <- removeBatchEffect(counts_b_CTF,design=treatment.design_b,covariates=batch.design_b)
write.table(counts_b_CTF_corr,"coex/blood/corrected_counts_CTF_blood.csv", sep = ",", quote=FALSE)

## DESeq2 VST
# Correct counts for variation not interesting
counts_b <- counts_b[,rownames(colData_b) ]
counts_b_corr <- removeBatchEffect(counts_b,design=treatment.design_b,covariates=batch.design_b)
write.table(counts_b_corr,"coex/blood/corrected_counts_blood.csv", sep = ",", quote=FALSE)

# PCA uncorrected
counts_b_log <- log2(counts_b+0.25)
counts_b_CTF_log <- log2(counts_b_CTF+0.25)

res.pca_b <- prcomp(t(counts_b_CTF_log), scale = F,center =T)

pdf(file = "coex/blood/pca_blood_CTF.pdf", width = 6, height = 4)
fviz_pca_ind(res.pca_b,
             geom = c("point"),
             col.ind = colData_b$Status, # color by groups
             addEllipses = T, # Concentration ellipses
             ellipse.level=0.95,
             legend.title = "Groups",
             repel = T
)
dev.off()

# PCA corrected
counts_b_corr_log <- log2(counts_b_corr+1)
counts_b_corr_log <- as.data.frame(counts_b_corr_log[complete.cases(counts_b_corr_log), ])
counts_b_CTF_corr_log <- log2(counts_b_CTF_corr+0.25)
counts_b_CTF_corr_log <- as.data.frame(counts_b_CTF_corr_log[complete.cases(counts_b_CTF_corr_log), ])

res.pca_b_c <- prcomp(t(counts_b_corr_log), scale = F,center =T)

pdf(file = "coex/blood/pca_corrected_blood_tissues.pdf", width = 6, height = 4)
fviz_pca_ind(res.pca_b_c,
             geom = c("point"),
             col.ind = colData_b$Tissue2, # color by groups
             addEllipses = T, # Concentration ellipses
             ellipse.level=0.95,
             legend.title = "Groups",
             repel = T
)
dev.off()

# Heatmap
pheatmap(cor(counts_b_corr_log), 
         annotation = colData_b[c("TissueBP","Status")],
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         show_colnames = T)

########################################################################################
### Co-expression GWENA / WGCNA ###
threads <- 7

counts_b_corr_f <- filter_low_var(t(counts_b_corr), pct = 0.7, type = "median")

netb <- build_net(counts_b_corr_f, cor_func = "spearman", n_threads = threads)

# CLR transformation of network
netb_clr <- netb
netb_clr$network <- clr(netb$network)

# Power selected :
netb$metadata$power

# Fit of the power law to data ($R^2$) :
fit_power_table <- netb$metadata$fit_power_table
fit_power_table[fit_power_table$Power == netb$metadata$power, "SFT.R.sq"]
plot(fit_power_table$Power, fit_power_table$SFT.R.sq)
plot(fit_power_table$Power, fit_power_table$mean.k.)

# without CLR
modulesb <- detect_modules(counts_b_corr_f, 
                          netb$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.9,
                          min_module_size = 30)
# with CLR
modulesb <- detect_modules(counts_b_corr_f, 
                           netb_clr$network, 
                           detailled_result = TRUE,
                           merge_threshold = 0.75,
                           min_module_size = 30)

# Number of modules pre/after merging:
length(unique(modulesb$modules_premerge))
length(unique(modulesb$modules))

# plot expression profiles of modules
genomic_idx <- match(rownames(colData_b[order(colData_b$Status),]), rownames(counts_b_corr_f))
expr_ordered  <- counts_b_corr_f[genomic_idx,]
plot_expression_profiles(expr_ordered, modulesb$modules$`1`)

# save modules in txt file
modules_df <- data.frame(Gene_ID=character(), Module=character(), stringsAsFactors=FALSE)
for (mod in names(modulesb$modules)) {
  x <- data.frame(Gene_ID=modulesb$modules[[mod]], stringsAsFactors=FALSE)
  x["Module"] <- rep(mod, length(modulesb$modules[[mod]]))
  modules_df <- rbind(modules_df, x)
}
write.table(modules_df,"coex/blood/modules_blood.txt", sep = "\t", quote=FALSE, row.names = FALSE)
table(modules_df$Module) # Number genes in modules

# plot merging of modules
layout_mod_merge <- plot_modules_merge(
  modules_premerge = modulesb$modules_premerge, 
  modules_merged = modulesb$modules)

#plot module sizes
pdf(file = "coex/blood/modules.pdf", width = 5, height = 3)
ggplot2::ggplot(data.frame(modulesb$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")
dev.off()

# Phenotype association
phenotype_association <- associate_phenotype(
  modulesb$modules_eigengenes, 
  colData_b[c("Status", "Condition")])

pdf(file = "coex/blood/phenotype_blood.pdf", width = 9, height = 4)
plot_modules_phenotype(phenotype_association, pvalue_th = 0.001)
dev.off()

write.table(phenotype_association$pval,"coex/blood/assoc_pval_blood.txt", sep = "\t", quote=FALSE)
write.table(phenotype_association$association,"coex/blood/assoc_r_blood.txt", sep = "\t", quote=FALSE)

# Module enrichments
enrichment <- bio_enrich(modulesb$modules[1:5], 
                         organism="oarambouillet",
                         correction_method = "fdr", 
                         domain_scope = "custom_annotated", 
                         custom_bg = colnames(counts_b_corr_f))

enrichment <- bio_enrich(modulesb$modules$`19`, 
                         organism="oarambouillet",
                         correction_method = "fdr", 
                         domain_scope = "custom_annotated", 
                         custom_bg = colnames(counts_b_corr_f),
                         significant = TRUE)

plot_enrichment(enrichment)

# Save enrichments to files
for (mod in names(modulesb$modules)) {
  gores <- c()
  enrichment <- bio_enrich(modulesb$modules[[mod]], 
                           organism="oarambouillet",
                           correction_method = "fdr", 
                           domain_scope = "custom_annotated", 
                           custom_bg = colnames(counts_b_corr_f),
                           significant = TRUE)
  
  if (length(rownames(enrichment$result)) > 1) {
    if (is.null(enrichment$result) == FALSE) {
      gores <- apply(enrichment$result,2,as.character)
      gores <- gores[, colnames(gores) != "parents"]
      gorestxt <- paste("coex/blood/gprofiler/ME", mod, ".txt", sep = "")
      write.table(gores, gorestxt, sep = "\t", quote=FALSE, row.names = FALSE)
    }
  }
}


# Hub genes
hubgenes <- get_hub_high_co(netb$network, modules = modulesb$modules, top_n = 10)
write.table(hubgenes,"coex/blood/hubmodules_blood.txt", sep = "\t", quote=FALSE, row.names = FALSE)

# Association false positive rate
# Correlations between random module eigengenes and treatment
pvals <- c()
for (i in 1:1000) {
  samp <- sample(nrow(t(counts_b_corr_f)),sample(30:1000, 1))
  sampe <- t(counts_b_corr_f)[samp,]
  eigengenes <- svd(sampe)$v
  t <- cor.test(eigengenes[,1], treatment.design_b[,"StatusTreatment"])
  pvals <- append(pvals, t$p.value)
}
# Select threshold p-value from plot
sum(pvals < 0.0001) / length(pvals)
fpr <- function(x) sum(pvals < x) / length(pvals)
p <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
fprs <- lapply(p, fpr)
plot(p, fprs, log="x")

########################################################################################
### Module preservation test ###
colData_b_ctr <- colData_b[colData_b["Status"]=="Control",]
colData_b_tre <- colData_b[colData_b["Status"]=="Treatment",]
cond_b_ctr <- counts_b_corr_f[row.names(colData_b_ctr),]
cond_b_tre <- counts_b_corr_f[row.names(colData_b_tre),]

expr_by_condb <- list("Control" = cond_b_ctr, "Treatment" = cond_b_tre)

net_by_condb <- lapply(expr_by_condb, build_net, cor_func = "spearman", 
                      n_threads = threads, keep_matrices = "both",
                      fit_cut_off = 0.85)

mod_by_condb <- mapply(detect_modules, expr_by_condb, 
                      lapply(net_by_condb, `[[`, "network"), 
                      MoreArgs = list(detailled_result = TRUE, merge_threshold = 0.85), 
                      SIMPLIFY = FALSE)

net_by_condb$Control$metadata$fit_power_table
net_by_condb$Control$metadata$power
net_by_condb$Treatment$metadata$fit_power_table
net_by_condb$Treatment$metadata$power

ggplot2::ggplot(data.frame(mod_by_condb$Control$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

#
comparisonb <- compare_conditions(expr_by_condb, 
                                 lapply(net_by_condb, `[[`, "adja_mat"), 
                                 lapply(net_by_condb, `[[`, "cor_mat"),  
                                 lapply(mod_by_condb, `[[`, "modules"), 
                                 pvalue_th = 0.01)

comparisonb$result$Control$Treatment$comparison
plot_comparison_stats(comparisonb$result$Control$Treatment$p.values)

########################################################################################
### Differential co-expression dcanr / z-score ###
# Prepare input metadata
dce_colData_b <- colData_b["Status"]
dce_colData_b["Condition"] <- 1
dce_colData_b$Condition[dce_colData_b$Status == "Treatment"] <- 2
dce_design_b <- dce_colData_b[["Condition"]]
names(dce_design_b) <- rownames(dce_colData_b)

# Calculate DC (select method and correlation metric)
z_scores <- dcScore(t(counts_b_corr_f), dce_design_b, dc.method = 'zscore', cor.method = 'spearman')
raw_p <- dcTest(z_scores, t(counts_b_corr_f), dce_design_b)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')

# Make DC network
dcnet <- dcNetwork(z_scores, adj_p, thresh = 5E-2) # filter by FDR
adjmat <- as_adj(dcnet, sparse = FALSE)
edgedf <- as_data_frame(dcnet, what = 'edges')

# Unique number of genes in DGC network
length(unique(c(edgedf[["from"]],edgedf[["to"]])))

# Save DGC network (for cytoscape)
write.table(edgedf,"coex/blood/network.txt", sep = "\t", quote=FALSE, row.names = FALSE)

# Sub-networks based on co-expression modules (for cytoscape)
edgemod = merge(edgedf, modules_df, by.x = "from", by.y = "Gene_ID")
edgemod = merge(edgemod, modules_df, by.x = "to", by.y = "Gene_ID")

for (mod in names(modulesb$modules)) {
  edgedfsub <- edgemod[edgemod["Module.x"] == mod | edgemod["Module.y"] == mod,]
  txt <- paste("coex/blood/subnet/network", mod, ".txt", sep = "")
  write.table(edgedfsub,txt, sep = "\t", quote=FALSE, row.names = FALSE)
}

# Transcription factors in subnetworks
tfs <- read.csv("references/TF_oaramb.csv")
tfs <- tfs[tfs$name!="None",]

edgedfsub <- edgemod[edgemod["Module.x"] == "16" | edgemod["Module.y"] == "16",]
subn_uniq <- unique(c(edgedfsub[["from"]],edgedfsub[["to"]]))
tfs[tfs$converted_alias %in% subn_uniq,]

for (mod in names(modulesb$modules)) {
  edgedfsub <- edgemod[edgemod["Module.x"] == mod | edgemod["Module.y"] == mod,]
  subn_uniq <- unique(c(edgedfsub[["from"]],edgedfsub[["to"]]))
  tfsmod <- tfs[tfs$converted_alias %in% subn_uniq,]
  txt <- paste("coex/blood/subnet/TF_network", mod, ".txt", sep = "")
  write.table(tfsmod,txt, sep = "\t", quote=FALSE, row.names = FALSE)
}

########################################################################################
save.image(file="coex/coex_atlas.RData")

########################################################################################
### LYMPH NODES ########################################################################
########################################################################################
load("coex/coex_atlas_lymph.RData")

## Test two normalisation methods: 1: Deseq2 VST 2: EdgeR CTF + asinh

# Design for covariate correction
design_t <- model.matrix(~Status + TissueBP, data=colData_t)
treatment.design_t <- design_t[,1:2]
batch.design_t <- design_t[,-(1:2)]

## EdgeR CTF + asinh
# CTF + asinh normalisation
rawcounts_t <- counts[,rownames(colData_t) ]
rawcounts_t <- rawcounts_t[rownames(counts_t), ]
norm_factors_t <- calcNormFactors(object = rawcounts_t, lib.size = colSums(rawcounts_t), method = "TMM")
counts_t_CTF <- sweep(rawcounts_t, 2, norm_factors_t, "/")
counts_t_CTF <- asinh(counts_t_CTF)
# Correct counts for variation not interesting
counts_t_CTF_corr <- removeBatchEffect(counts_t_CTF,design=treatment.design_t,covariates=batch.design_t)
write.table(counts_t_CTF_corr,"coex/lymph/corrected_counts_CTF_lymph.csv", sep = ",", quote=FALSE)

## DESeq2 VST
# Correct counts for variation not interesting
counts_t <- counts_t[,rownames(colData_t) ]
counts_t_corr <- removeBatchEffect(counts_t,design=treatment.design_t,covariates=batch.design_t)
write.table(counts_t_corr,"coex/lymph/corrected_counts_lymph.csv", sep = ",", quote=FALSE)

# PCA uncorrected
counts_t_log <- log2(counts_t+0.25)
counts_t_CTF_log <- log2(counts_t_CTF+0.25)

res.pca_t <- prcomp(t(counts_t_log), scale = F,center =T)

pdf(file = "coex/lymph/pca_lymph_bp.pdf", width = 6, height = 4)
fviz_pca_ind(res.pca_t,
             geom = c("point"),
             col.ind = colData_t$BioProject, # color by groups
             addEllipses = T, # Concentration ellipses
             ellipse.level=0.95,
             legend.title = "Groups",
             repel = T
)
dev.off()

# PCA corrected
counts_t_corr_log <- log2(counts_t_corr+1)
counts_t_corr_log <- as.data.frame(counts_t_corr_log[complete.cases(counts_t_corr_log), ])
counts_t_CTF_corr_log <- log2(counts_t_CTF_corr+0.25)
counts_t_CTF_corr_log <- as.data.frame(counts_t_CTF_corr_log[complete.cases(counts_t_CTF_corr_log), ])

res.pca_t_c <- prcomp(t(counts_t_corr_log), scale = F,center =T)

pdf(file = "coex/lymph/pca_corrected_lymph_bp.pdf", width = 6, height = 4)
fviz_pca_ind(res.pca_t_c,
             geom = c("point"),
             col.ind = colData_t$Status, # color by groups
             addEllipses = T, # Concentration ellipses
             ellipse.level=0.95,
             legend.title = "Groups",
             repel = T
)
dev.off()

# Heatmap
pheatmap(cor(counts_t_corr_log), 
         annotation = colData_t[c("BioProject","Status")],
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         show_colnames = T)

########################################################################################
### Co-expression GWENA / WGCNA ###
threads <- 7

counts_t_corr_f <- filter_low_var(t(counts_t_corr), pct = 0.7, type = "median")

nett <- build_net(counts_t_corr_f, fit_cut_off = 0.85, cor_func = "spearman", n_threads = threads)

# CLR transformation of network
nett_clr <- nett
nett_clr$network <- clr(nett$network)

# Power selected :
nett$metadata$power

# Fit of the power law to data ($R^2$)
fit_power_table <- nett$metadata$fit_power_table
fit_power_table[fit_power_table$Power == nett$metadata$power, "SFT.R.sq"]
plot(fit_power_table$Power, fit_power_table$SFT.R.sq)
plot(fit_power_table$Power, fit_power_table$mean.k.)

# without CLR
modulest <- detect_modules(counts_t_corr_f, 
                           nett$network, 
                           detailled_result = TRUE,
                           merge_threshold = 0.85,
                           min_module_size = 30)

# with CLR
modulest <- detect_modules(counts_t_corr_f, 
                           nett_clr$network, 
                           detailled_result = TRUE,
                           merge_threshold = 0.75,
                           min_module_size = 30)

# Number of modules pre/after merging: 
length(unique(modulest$modules_premerge))
length(unique(modulest$modules))

# plot expression profiles of modules
genomic_idx <- match(rownames(colData_t[order(colData_t$Status),]), rownames(counts_t_corr_f))
expr_ordered  <- counts_t_corr_f[genomic_idx,]
plot_expression_profiles(expr_ordered, modulest$modules$`1`)

# save modules in txt file
modules_df <- data.frame(Gene_ID=character(), Module=character(), stringsAsFactors=FALSE)
for (mod in names(modulest$modules)) {
  x <- data.frame(Gene_ID=modulest$modules[[mod]], stringsAsFactors=FALSE)
  x["Module"] <- rep(mod, length(modulest$modules[[mod]]))
  modules_df <- rbind(modules_df, x)
}
write.table(modules_df,"coex/lymph/modules_lymph.txt", sep = "\t", quote=FALSE, row.names = FALSE)
table(modules_df$Module)

# plot merging of modules
layout_mod_merge <- plot_modules_merge(
  modules_premerge = modulest$modules_premerge, 
  modules_merged = modulest$modules)

#plot module sizes
pdf(file = "coex/lymph/modules.pdf", width = 5, height = 3)
ggplot2::ggplot(data.frame(modulest$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")
dev.off()

# Module enrichment
enrichment <- bio_enrich(modulest$modules[13:20], 
                         organism="oarambouillet",
                         correction_method = "fdr", 
                         domain_scope = "custom_annotated", 
                         custom_bg = colnames(counts_t_corr_f))

enrichment <- bio_enrich(modulest$modules$`31`, 
                         organism="oarambouillet",
                         correction_method = "fdr", 
                         domain_scope = "custom_annotated", 
                         custom_bg = colnames(counts_t_corr_f),
                         significant = TRUE)

plot_enrichment(enrichment)

# Save enrichments to files
for (mod in names(modulest$modules)) {
  gores <- c()
  enrichment <- bio_enrich(modulest$modules[[mod]], 
                           organism="oarambouillet",
                           correction_method = "fdr", 
                           domain_scope = "custom_annotated", 
                           custom_bg = colnames(counts_t_corr_f),
                           significant = TRUE)
  
  if (length(rownames(enrichment$result)) > 1) {
    if (is.null(enrichment$result) == FALSE) {
      gores <- apply(enrichment$result,2,as.character)
      gores <- gores[, colnames(gores) != "parents"]
      gorestxt <- paste("coex/lymph/gprofiler/ME", mod, ".txt", sep = "")
      write.table(gores, gorestxt, sep = "\t", quote=FALSE, row.names = FALSE)
    }
  }
}


# Phenotype association
phenotype_association <- associate_phenotype(
  modulest$modules_eigengenes, 
  colData_t[c("Status", "Condition")])

pdf(file = "coex/lymph/phenotype_lymph.pdf", width = 9, height = 3)
plot_modules_phenotype(phenotype_association, pvalue_th = 0.001)
dev.off()

write.table(phenotype_association$pval,"coex/lymph/assoc_pval_lymph.txt", sep = "\t", quote=FALSE)
write.table(phenotype_association$association,"coex/lymph/assoc_r_lymph.txt", sep = "\t", quote=FALSE)

# Association false positive rate
# Correlations between random module eigengenes and treatment
pvals <- c()
for (i in 1:1000) {
  samp <- sample(nrow(t(counts_t_corr_f)),sample(30:1000, 1))
  sampe <- t(counts_t_corr_f)[samp,]
  eigengenes <- svd(sampe)$v
  t <- cor.test(eigengenes[,1], treatment.design_t[,"StatusTreatment"])
  pvals <- append(pvals, t$p.value)
}
# Select threshold p-value from plot
sum(pvals < 0.0001) / length(pvals)
fpr <- function(x) sum(pvals < x) / length(pvals)
p <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
fprs <- lapply(p, fpr)
plot(p, fprs, log="x")

# Hub genes
hubgenes <- get_hub_high_co(nett$network, modules = modulest$modules, top_n = 10)
hubgenes

########################################################################################
# Module preservation test
colData_t_ctr <- colData_t[colData_t["Status"]=="Control",]
colData_t_tre <- colData_t[colData_t["Status"]=="Treatment",]
cond_b_ctr <- counts_t_corr_f[row.names(colData_t_ctr),]
cond_b_tre <- counts_t_corr_f[row.names(colData_t_tre),]

expr_by_cond <- list("Control" = cond_b_ctr, "Treatment" = cond_b_tre)

net_by_cond <- lapply(expr_by_cond, build_net, cor_func = "spearman", 
                      n_threads = threads, keep_matrices = "both",
                      fit_cut_off = 0.85)

mod_by_cond <- mapply(detect_modules, expr_by_cond, 
                      lapply(net_by_cond, `[[`, "network"), 
                      MoreArgs = list(detailled_result = TRUE, merge_threshold = 0.8), 
                      SIMPLIFY = FALSE)

net_by_cond$Control$metadata$fit_power_table
net_by_cond$Control$metadata$power
net_by_cond$Treatment$metadata$fit_power_table
net_by_cond$Treatment$metadata$power

ggplot2::ggplot(data.frame(mod_by_cond$Control$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

#
comparison <- compare_conditions(expr_by_cond, 
                                 lapply(net_by_cond, `[[`, "adja_mat"), 
                                 lapply(net_by_cond, `[[`, "cor_mat"),  
                                 lapply(mod_by_cond, `[[`, "modules"), 
                                 pvalue_th = 0.05)

comparison$result$Control$Treatment$comparison
plot_comparison_stats(comparison$result$Control$Treatment$p.values)

########################################################################################
### Differential co-expression dcanr / z-score ###
# Prepare input metadat
dce_colData_t <- colData_t["Status"]
dce_colData_t["Condition"] <- 1
dce_colData_t$Condition[dce_colData_t$Status == "Treatment"] <- 2
dce_design_t <- dce_colData_t[["Condition"]]
names(dce_design_t) <- rownames(dce_colData_t)
table(dce_design_t)

# Calculate DC (select method and correlation metric)
z_scores <- dcScore(t(counts_t_corr_f), dce_design_t, dc.method = 'zscore', cor.method = 'spearman')
raw_p <- dcTest(z_scores, t(counts_t_corr_f), dce_design_t)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')

# Make DC network
dcnet <- dcNetwork(z_scores, adj_p, thresh = 1E-3) # filter by FDR 1E-3
adjmat <- as_adj(dcnet, sparse = FALSE)
edgedf <- as_data_frame(dcnet, what = 'edges')

# Unique number of genes in DGC network
length(unique(c(edgedf[["from"]],edgedf[["to"]])))

# Save DGC network
write.table(edgedf,"coex/lymph/network_lymph.txt", sep = "\t", quote=FALSE, row.names = FALSE)

# Sub-networks based on co-expression modules
edgemod = merge(edgedf, modules_df, by.x = "from", by.y = "Gene_ID")
edgemod = merge(edgemod, modules_df, by.x = "to", by.y = "Gene_ID")

for (mod in names(modulest$modules)) {
  edgedfsub <- edgemod[edgemod["Module.x"] == mod | edgemod["Module.y"] == mod,]
  txt <- paste("coex/lymph/subnet/network", mod, ".txt", sep = "")
  write.table(edgedfsub,txt, sep = "\t", quote=FALSE, row.names = FALSE)
}

# Transcription factors in subnetworks
tfs <- read.csv("references/TF_oaramb.csv")
tfs <- tfs[tfs$name!="None",]

for (mod in names(modulest$modules)) {
  edgedfsub <- edgemod[edgemod["Module.x"] == mod | edgemod["Module.y"] == mod,]
  subn_uniq <- unique(c(edgedfsub[["from"]],edgedfsub[["to"]]))
  tfsmod <- tfs[tfs$converted_alias %in% subn_uniq,]
  txt <- paste("coex/lymph/subnet/TF_network", mod, ".txt", sep = "")
  write.table(tfsmod,txt, sep = "\t", quote=FALSE, row.names = FALSE)
}
############################################################################
save.image(file="coex/coex_atlas_lymph.RData")


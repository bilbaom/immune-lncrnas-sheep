#!/usr/bin/R

#============================================== CAGE data analysis using CAGEfightR pipeline ==============================================
# The bigWig files then were imported to CAGEfightR within R environment
# The sheep genome BSgenome object was created following instructions of BSgenome package for costum genomes:
# https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
# https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf

require(CAGEfightR)
require(magrittr)
require(GenomicRanges)
require(BSgenome)
require(MultiAssayExperiment)
require(SummarizedExperiment)
require(Biostrings)
require(GenomicFeatures)
require(BiocParallel)
require(InteractionSet)
require(Gviz)
require(BSgenome.Oaries.NCBI.Ramb1.0)
require(tidyverse)
require(FactoMineR)
require(factoextra)
require(reshape2)
require(circlize)
require(RColorBrewer)
require(randomcoloR) #!!
require(philentropy)
require(bioDist)
require(patchwork)
require(tidyverse)
require(scales)
require(tximport)
require(Repitools)
require(Hmisc)
require(viridis)
require(lme4)
require(lmerTest)
setwd("/home/labo/datuak/atlas/cage")

#Reading the bw files in to CAGEfightR

#Using the name of the BAM files to reorder the tissues tally numbers
inputfiles<-list.files(path = "bam",pattern = ".bam$",full.names = TRUE,include.dirs =TRUE)
inputnames<-gsub(x = str_split_fixed(inputfiles,"/",2)[,2],pattern = ".bam",replacement = "")

tally <- data.frame (Sample  = inputnames,
                     Tissue = c("alveolar-macrophage", "tonsil-palatine", "spleen","lymph node-mesenteric","lymph node-prescapular"))

# The BSgenome object for both assemblies were tested an imported as separate object in R. 
Ramb1_Seqinfo<-Seqinfo(seqnames = BSgenome.Oaries.NCBI.Ramb1.0@seqinfo@seqnames,
        seqlengths = BSgenome.Oaries.NCBI.Ramb1.0@seqinfo@seqlengths,
        isCircular = BSgenome.Oaries.NCBI.Ramb1.0@seqinfo@is_circular,
        genome = BSgenome.Oaries.NCBI.Ramb1.0@seqinfo@genome)

#Reading BigWig files into CAGEfightR
bw_plus <-
  dir(
    path = "bam",
    pattern = ".plus.bw",
    all.files = T,
    full.names = T,
    include.dirs = T
  )
bw_minus <-
  dir(
    path = "bam",
    pattern = ".minus.bw",
    all.files = T,
    full.names = T,
    include.dirs = T
  )

#Confirming equal number of plus and minus tissues within the folder
bw_minus<-bw_minus[match(str_replace_all(bw_plus,pattern = ".plus.bw",replacement = ""),str_replace_all(bw_minus,pattern=".minus.bw",""))]

#Creating the BigWig list objects
bw_plus<-BigWigFileList(bw_plus)
bw_minus<-BigWigFileList(bw_minus)
names(bw_minus) <- inputnames
names(bw_plus) <- inputnames

#As the conversion of the bigWigs are done using 1 indexed coordinates the legth of the BSgenome is 1 bp short of the bw files. Fixed it manually.

BSgenome.Oaries.NCBI.Ramb1.0@seqinfo@seqlengths<-as.integer(Ramb1_Seqinfo@seqlengths+1)

CTSSs<-quantifyCTSSs(plusStrand = bw_plus,
                     minusStrand = bw_minus,
                     genome = SeqinfoForBSGenome(BSgenome.Oaries.NCBI.Ramb1.0))



CTSSs<-CTSSs %>% calcTPM() %>% calcPooled()
TCs<-quickTSSs(CTSSs)

TSSs <- TCs %>%
  calcTPM() %>%
  subsetBySupport(inputAssay="counts",
                  unexpressed=10, #minimum TPM requirement
                  minSamples=0) %>%  #Keep all 
  calcTPM()

# Same issue with the chromosome lengths of the txdb object
chr_seqlength<-read.csv("bam/chromsizes.txt",header = F,stringsAsFactors = F,sep = ' ')
#chr_seqlength$V2<-chr_seqlength$V2+1

#Making the Seqinfo object
txdb_fix<-Seqinfo(seqnames = as.character(chr_seqlength$V1),
                  seqlengths=as.integer(chr_seqlength$V2),
                  isCircular = ifelse(chr_seqlength$V1=="NC_001941.1",TRUE,FALSE))
				  
#Creating sheep TxDB
#### CHANGE CHROMOSOME NAMES FROM ENSEMBL TO NCBI!!!!!
txdb<-makeTxDbFromGFF(file = "/home/labo/datuak/atlas/cage/cagefight.gtf",
                      format = "auto",
                      dataSource = NA,
                      organism = "Ovis aries",
                      chrominfo = txdb_fix, #forcing the seqlength info by chrominfo
                      taxonomyId = 9940,
                      dbxrefTag = TRUE)

#Matching the length of chromosomes
seqlevels(txdb)<-seqlevels(TCs)

TSSs <- assignTxType(TSSs,txModels = txdb,swap="thick", tssUpstream = 500, tssDownstream = 500)
TSSs<-assignTxID(object = TSSs,txModels = txdb,swap="thick")
TSSs <- assignGeneID(TSSs, txdb)
TSSs <- assignMissingID(TSSs, outputColumn = "geneID", prefix = "NOVELG")
TSSs <- assignMissingID(TSSs, outputColumn = "txID", prefix = "NOVELT")

# Results to dataframe
TSSs_df<-TSSs %>% rowRanges() %>% data.frame()
write.table(TSSs_df,"tss_cage.csv", sep = ",", quote=FALSE)

#Bidirectional clustering
BCs<-quickEnhancers(CTSSs)
BCs<-subsetBySupport(BCs, inputAssay="counts", unexpressed=10, minSamples=0)
BCs<-assignTxType(BCs, txModels = txdb, swap="thick", tssUpstream = 500, tssDownstream = 500)
BCs<-assignTxID(BCs, txModels = txdb, swap="thick")
BCs <- assignGeneID(BCs, txdb)
BCs <- assignMissingID(BCs, outputColumn = "geneID", prefix = "NOVELG")
BCs <- assignMissingID(BCs, outputColumn = "txID", prefix = "NOVELT")

# Results to dataframe
BCs_df<-BCs %>% rowRanges() %>% data.frame()
write.table(BCs_df,"bidirectional_cage.csv", sep = ",", quote=FALSE)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

Enhancers <- BCs %>%
  calcTPM() %>%
  subsetBySupport(inputAssay="counts",
                  unexpressed=10, #minimum TPM requirement
                  minSamples=0) %>% #The 2/3rd representation criteria. 
  calcTPM()

Enhancers_df<-Enhancers %>% rowRanges() %>% data.frame()


BCs_full<-subsetBySupport(BCs, inputAssay="counts", unexpressed=10, minSamples=1)

##saving point
con <- pipe("xz -T8 -6 -e > crf_BCs_object.xz", "wb")
save(BCs, file = con); close(con)
rm(con)

con <- pipe("xz -T8 -6 -e > crf_Enhancers_object.xz", "wb")
save(Enhancers, file = con); close(con)
rm(con)

####################################################################################################
#Calculating co-expression correlation coefficience
TSSs$totalTags <- NULL
BCs$totalTags <- NULL
 
CAGEclusters <- combineClusters(object1 = TSSs,
                                 object2 = BCs)
 
 rowRanges(CAGEclusters)$clusterType <-
   factor(ifelse(strand(CAGEclusters) == "*",
                 "enhancer", "TSS"),
				 levels = c("TSS", "enhancer"))
 links <- findLinks(
   CAGEclusters,
   inputAssay = "counts",
   directional = "clusterType",
   method = "kendall"
) %>% subset(p.value < 0.05)

#Storing various input dataframe for RNA-Seq overlay workup. 

links_df<- regions(links) %>% as.data.frame()
TSSs_df<-TSSs %>% rowRanges() %>% data.frame()
BCs_df<-BCs %>% rowRanges() %>% data.frame()
TSSs_pca_input<-TSSs@assays@.xData$data$TPM
TSSs_full_pca_input<-TSSs_full@assays@.xData$data$TPM
BCs_pca_input<-Enhancers@assays@.xData$data$TPM



#Tissue specific work up 


TSS_PCA<-PCA(TSSs_pca_input,scale.unit =T,ncp = 5,graph = F)
BC_PCA<-PCA(BCs_pca_input,scale.unit = T,ncp = 5,graph = F)


cage_write <- function(df,tissue) {
  tmp <- data.frame(
    cluster = row.names(TSSs_pca_input),
    tissue = TSSs_pca_input[,tissue],
    row.names = 1:nrow(TSSs_pca_input)
  ) %>%
    filter(tissue != 0)
  tmp1 <- tmp$cluster %>% str_split_fixed(":", 2)
  tmp2 <- tmp1[, 2] %>% str_split_fixed(";", 2)
  tmp3 <- tmp2[, 1] %>% str_split_fixed("-", 2)
  ready <-
    cbind.data.frame(
      chr = tmp1[, 1],
      start = tmp3[, 1],
      end = tmp3[, 2],
      strand = tmp2[, 2],
      TPM = tmp$tissue
    )
  write.table(
    ready,
    paste0("output_CAGEfightR/", tissue, "_TSSs_tpm.tsv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

for (i in seq_along(colnames(TSSs_pca_input))){
  cage_write(TSSs_pca_input,colnames(TSSs_pca_input)[i])
}
  

cage_write2 <- function(df,tissue) {
  tmp <- data.frame(
    cluster = row.names(BCs_pca_input),
    tissue = BCs_pca_input[,tissue],
    row.names = 1:nrow(BCs_pca_input)
  ) %>%
    filter(tissue != 0)
  tmp1 <- tmp$cluster %>% str_split_fixed(":", 2)
  tmp2 <- tmp1[, 2] %>% str_split_fixed("-", 2)
  ready <-
    cbind.data.frame(
      chr = tmp1[, 1],
      start = tmp2[, 1],
      end = tmp2[, 2],
      #strand = tmp2[, 2],
      TPM = tmp$tissue
    )
  write.table(
    ready,
    paste0("output_CAGEfightR/", tissue, "_BCs_tpm.tsv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

for (i in seq_along(colnames(BCs_pca_input))){
  cage_write2(BCs_pca_input,colnames(BCs_pca_input)[i])
}

#Calulating MI distance between tissues given BC or TSS cluster expression information. 
BCs_invert<-t(BCs_pca_input)
BCs_promoters_invert<-t(Enhancers %>% subset(txType=="promoter") %>% assays() %>% .$TPM)
BCs_CDS_invert<-t(Enhancers %>% subset(txType=="CDS") %>% assays() %>% .$TPM)

TSSs_invert<-t(TSSs_pca_input)
TSSs_promoters_invert<-t(TSSs %>% subset(txType=="promoter") %>% assays() %>% .$TPM)
TSSs_CDS_invert<-t(TSSs %>% subset(txType=="CDS") %>% assays() %>% .$TPM)


TSSs_MI<-MIdist(TSSs_invert)
TSSs_MI_promoters<-MIdist(TSSs_promoters_invert)
TSSs_MI_CDS<-MIdist(TSSs_CDS_invert)


BCs_MI<-MIdist(BCs_invert)
BCs_MI_promoters<-MIdist(BCs_promoters_invert)
BCs_MI_CDS<-MIdist(BCs_CDS_invert)

# The annotation metrics of tissue specific CAGE tag clusters
#TSS set
seqlevels(txdb)<-seqlevels(TSSs_full)
TSSs_full <- assignTxType(TSSs_full,
                          txModels = txdb,
                          swap="thick")
TSSs_full<-assignTxID(object = TSSs_full,
                      txModels = txdb,
                      swap="thick")
TSSs_full_promoters<-TSSs_full %>% subset(txType=="promoter") %>% subset(support==1)
TSSs_full_promoters_tpm<-assay(TSSs_full_promoters,"TPM")
TSSs_full_promoters_df<-rowRanges(TSSs_full_promoters) %>% as.data.frame()
rownames(TSSs_full_promoters_df)<-1:nrow(TSSs_full_promoters_df)
TSSs_full_promoters_tpm<-data.frame(thick.names=rownames(TSSs_full_promoters_tpm),TSSs_full_promoters_tpm)
rownames(TSSs_full_promoters_tpm)<-1:nrow(TSSs_full_promoters_tpm)

specific_metrics<-data.frame(tissue=NULL,total=NULL,TPMabove1=NULL)
for (col in seq_along(TSSs_full_promoters_tpm)[-1]){
  tpm=inner_join(TSSs_full_promoters_tpm[which(TSSs_full_promoters_tpm[col]!=0),c(1,col)],TSSs_full_promoters_df[c(10,14)],"thick.names") 
  above1=tpm %>% filter(tpm[2]>1) %>% count(geneID)
  specific=cbind(tissue=colnames(TSSs_full_promoters_tpm)[col],total=length(tpm$geneID),TPMabove1=length(above1$geneID))
  specific_metrics<-rbind.data.frame(specific_metrics,specific)
}
specific_metrics[-1]<-map_df(specific_metrics[-1], function(x) x<-as.numeric(as.character(x)))

#BC set
BCs_full<-BCs_full %>% calcTPM()
BCs_full_promoters<-BCs_full %>% subset(txType=="promoter") %>% subset(support==1)

BCs_full_promoters_tpm<-assay(BCs_full_promoters,"TPM")
BCs_full_promoters_df<-rowRanges(BCs_full_promoters) %>% as.data.frame()
rownames(BCs_full_promoters_df)<-1:nrow(BCs_full_promoters_df)
BCs_full_promoters_tpm<-data.frame(thick.names=rownames(BCs_full_promoters_tpm),BCs_full_promoters_tpm)
rownames(BCs_full_promoters_tpm)<-1:nrow(BCs_full_promoters_tpm)

specific_metrics2<-data.frame(tissue=NULL,total=NULL,TPMabove1=NULL)
for (col in seq_along(BCs_full_promoters_tpm)[-1]){
  tpm=inner_join(BCs_full_promoters_tpm[which(BCs_full_promoters_tpm[col]!=0),c(1,col)],BCs_full_promoters_df[c(10,16)],"thick.names") 
  above1=tpm %>% filter(tpm[2]>1) %>% count(geneID)
  specific=cbind(tissue=colnames(BCs_full_promoters_tpm)[col],total=length(tpm$geneID),TPMabove1=length(above1$geneID))
  specific_metrics2<-rbind.data.frame(specific_metrics2,specific)
}
specific_metrics2[-1]<-map_df(specific_metrics2[-1], function(x) x<-as.numeric(as.character(x)))


#==============================================RNA-Seq overlay==============================================

dir_cage_tpm<-dir(path = "~/output_CAGEfightR",pattern = ".tsv",all.files = T,full.names = T,include.dirs = T)
dir_cage_tpm_tss<-dir_cage_tpm[str_detect(string = dir_cage_tpm,pattern = "TSSs")]
dir_cage_tpm_bc<-dir_cage_tpm[str_detect(string = dir_cage_tpm,pattern = "BC")]

#Reading the CAGE CTPMs in 
cage_tpm_tss <- list()
for (i in seq_along(dir_cage_tpm_tss)){
        name<-basename(dir_cage_tpm_tss[i]) %>% str_replace_all(pattern = "_TSSs_tpm.tsv","")
        tmp<-read_tsv(dir_cage_tpm_tss[i]) %>% janitor::clean_names()
        cage_tpm_tss[[i]]<-tmp
        names(cage_tpm_tss)[i]<-name
}

#Reading the RNA-Seq kallisto TPMs in

tsv_list<-dir(path = "~/input",pattern =".tsv",all.files = T,full.names = T,recursive = T,include.dirs = T)
tsv_list_OAR<-subset(tsv_list,subset = grepl(pattern = "OAR",x = tsv_list))
tsv_list_RAMB<-subset(tsv_list,subset = grepl(pattern = "RAMB",x = tsv_list))

input_OAR<-list()
for (i in seq_along(tsv_list_OAR)){
   name<-dirname(path = tsv_list_OAR[i]) %>% str_split_fixed(pattern = "/",n=8) %>% .[7]
   input_OAR[[i]]<-map_df(.x = tsv_list_OAR[i],.f = function(x) read.csv(x,header = T,sep = "\t",stringsAsFactors = F)) %>% filter(tpm>1)
   names(tsv_list_OAR)[i]<-names(input_OAR)[[i]]<-name
   rm(name,i)
 }
input_RAMB<-list()
for (i in seq_along(tsv_list_RAMB)){
  name<-dirname(path = tsv_list_RAMB[i]) %>% str_split_fixed(pattern = "/",n=8) %>% .[7]
  input_RAMB[[i]]<-map_df(.x = tsv_list_RAMB[i],.f = function(x) read.csv(x,header = T,sep = "\t",stringsAsFactors = F))
  names(tsv_list_RAMB)[i]<-names(input_RAMB)[[i]]<-name
  rm(name,i)
}

#transcript gene name map was produed from the GFF file. 
RAMB.tx2gene <-
  read.csv(
    "~/tx2gene.RAMB",
    sep = " ",
    header = F,
    stringsAsFactors = F
  )

tx2gene_RAMB <-RAMB.tx2gene %>%
  rename(GENEID=V2,TXID=V1) %>%
  select(TXID, GENEID) %>%
  as_tibble()

tx2gene_RAMB<-tx2gene_RAMB[!duplicated(tx2gene_RAMB$TXID),]

import_OAR <-
   tximport(
     files = tsv_list_OAR,
     type = "kallisto",
     tx2gene = tx2gene_OAR,
     ignoreTxVersion = T,
     ignoreAfterBar = F
   )

import_RAMB <-
  tximport(
    files = tsv_list_RAMB,
    type = "kallisto",
    tx2gene = tx2gene_RAMB,
    txOut = T,
    ignoreTxVersion = T,
    ignoreAfterBar = T
  )

RNA<-cbind.data.frame(txID=rownames(import_RAMB$abundance),
            import_RAMB$abundance) %>%
  separate(col = txID,into = c("chr","txID"),sep = ".1_",remove=T,extra = "merge") %>%
  separate(col = txID,into = c("type","txID"),sep = "_",remove=T,extra = "merge") %>%
  mutate(txID=str_remove(txID,pattern = "_\\d+$")) %>%
  dplyr::select(-chr,-type)

tally_RNA<-read_csv("~/tally_RNA-Seq.csv")
new<-unlist(tally_RNA[match(table = tally_RNA$run_accession,x = colnames(RNA)[-1]),2])
colnames(RNA)[-1]<-new

CAGE_short_list <- list()
CAGE_long_list <- list()
CAGE_novel_list <- list()
CAGE_annot_list <- list()
for (i in seq_along(cage_tpm_tss)) {
  tmp <- GenomicRanges::GRanges(
    seqnames = cage_tpm_tss[[i]]$chr,
    ranges = IRanges::IRanges(start = cage_tpm_tss[[i]]$start,
                              end = cage_tpm_tss[[i]]$end),
    strand = cage_tpm_tss[[i]]$strand,
    score = cage_tpm_tss[[i]]$tpm
  )
  CAGE_short_list[i] <- mergeByOverlaps(
    subject = tmp,
    #query = transcripts(txdb,c("tx_name","GENEID")),
    query = promoters(
      txdb,
      upstream = 25,
      downstream = 25,
      use.names = T,
      c("tx_name", "GENEID")
    ),
    maxgap = 25,
    type = "any"
  )
  CAGE_long_list[i] <- mergeByOverlaps(
    subject = tmp,
    query = promoters(
      txdb,
      upstream = 25,
      downstream = 25,
      use.names = T,
      c("tx_name", "GENEID")
    ),
    maxgap = 375,
    type = "any"
  )
  
  combo <-
    GRangesList(
      "short50bp" = CAGE_short_list[[i]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `,
      "long400bp" = CAGE_long_list[[i]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `)
  
  CAGE_novel_list[[i]] <- subsetByOverlaps(x = tmp,ranges = unlist(combo),invert = TRUE)
  CAGE_annot_list[[i]]<-subsetByOverlaps(x = tmp,ranges = CAGE_short_list[[i]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `)
  
  names(CAGE_short_list)[i] <- names(cage_tpm_tss[i])
  names(CAGE_long_list)[i] <- names(cage_tpm_tss[i])
  names(CAGE_novel_list)[i] <- names(cage_tpm_tss[i])
  names(CAGE_annot_list)[i]<-names(cage_tpm_tss[i])
  
  rm(tmp)
}

#Selecting the highest CTPM tag as the corresponding to the RNA-Seq transcript TPM. 

CAGE<-list()
for (i in seq_along(CAGE_short_list)) {
  tmp <-
    cbind.data.frame(
      cage = CAGE_short_list[[i]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `,
      tpm = CAGE_short_list[[i]]$score
    ) %>%
    select(cage.tx_name, tpm)%>%
      group_by(cage.tx_name) %>%
      filter(tpm == max(tpm)) #Only CAGE tag with the highest expression is selected
  colnames(tmp)[2] <- names(CAGE_short_list[i])
  CAGE[[i]]<-tmp %>% unique.data.frame() 
  names(CAGE)[i]<-names(cage_tpm_tss[i])
}

#Duplicated transcript names had to be removed from the CAGE input for it to pass through data.table merge operation. 
tx_reg<-data.frame(cage.tx_name=transcripts(txdb)$tx_name)
for (i in seq_along(CAGE)) {
  require(data.table)
  lhs<-as.data.table(tx_reg)
  rhs<-as.data.table(CAGE[[i]])
  tx_reg <- merge(lhs,rhs)
}
colnames(tx_reg)[1]<-"txID"

rna_reg<-subset(RNA,select=which(!duplicated(names(RNA))))
  
rna_cage_inner<-inner_join(
  reshape2::melt(rna_reg, value.name = "rna_tpm", id.vars = "txID"),
  reshape2::melt(
    tx_reg %>% drop_na(),
    value.name = "cage_tpm",
    id.vars = "txID"
  ),
  by = c("txID", "variable")
) %>% na.omit()


##Modelling the correlation between RNA-Seq and CAGE expression TPM in all transcripts

cage_rna_fit<-lmer(
  formula = log10(cage_tpm + 1) ~ log10(rna_tpm + 1) + (log10(rna_tpm + 1)|variable),
  rna_cage_inner,
  REML = F,
  control =  lmerControl(optimizer = "bobyqa")
)

summary(cage_rna_fit)

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log10(cage_tpm + 1) ~ log10(rna_tpm + 1) + (log10(rna_tpm + 1) |      variable)
#    Data: rna_cage_inner
# Control: lmerControl(optimizer = "bobyqa")
# 
#       AIC       BIC    logLik  deviance  df.resid 
#  447436.6  447500.2 -223712.3  447424.6    298058 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.1011 -0.6709 -0.0500  0.6161  5.7603 
# 
# Random effects:
#  Groups   Name               Variance Std.Dev. Corr 
#  variable (Intercept)        0.031965 0.17879       
#           log10(rna_tpm + 1) 0.004859 0.06971  -0.63
#  Residual                    0.262226 0.51208       
# Number of obs: 298064, groups:  variable, 52
# 
# Fixed effects:
#                    Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)         1.55202    0.02483 52.00034    62.5   <2e-16 ***
# log10(rna_tpm + 1)  0.18429    0.00975 51.70708    18.9   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#             (Intr)
# lg10(rn_+1) -0.630

novels<-data.frame(tissue=NA,novel_pct=NA,short=NA,long=NA,total=NA)
noveler<-function(TISSUE){
  ratio<-length(CAGE_novel_list[[TISSUE]]$score)/length(cage_tpm_tss[[TISSUE]]$tpm)
  tmp<-cbind(tissue=TISSUE,
             novel_pct=ratio,
             short=CAGE_short_list[[TISSUE]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `$GENEID %>% unlist %>% unique() %>% length(),
             long=CAGE_long_list[[TISSUE]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `$GENEID %>% unlist %>% unique() %>% length(),
             total=length(cage_tpm_tss[[TISSUE]]$tpm))
  novels<<-rbind.data.frame(tmp,novels,stringsAsFactors = F) %>% na.omit()
}
for (p in seq_along(tissue_reg)){
  noveler(tissue_reg[p])
}

novels[-1]<-lapply(novels[-1],as.numeric)
write_tsv(path="output/txdb_metrics_CAGE.tsv",novels)

#============================================== WGBS overlay ==============================================
#Processing WGBS data from RUMA and MUSLD

CAGE_writer<-function(LIST,IN,SET){
  NAME<-names(LIST)[IN]
  OUT<-data.frame(chr=seqnames(LIST[[IN]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `),
             start=start(LIST[[IN]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `),
             end=end(LIST[[IN]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `),
             width=width(LIST[[IN]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `),
             score=LIST[[IN]]$score,
             strand=strand(LIST[[IN]]$`promoters(txdb, upstream = 25, downstream = 25, use.names = T, `))
  assign(x = paste0(NAME),value = OUT)
  write_delim(x = get(NAME),
              path= paste0("output/",SET,"_CAGE/CAGE_",SET,"_",NAME,".bed"),
              delim = '\t')
}

novel_writer<-function(LIST,IN){
  NAME<-names(LIST)[IN]
  OUT<-annoGR2DF(LIST[[IN]])
  colnames(OUT)[1]<-"#chr"
  assign(x=NAME,value = OUT)
  write_tsv(x=get(NAME),path = paste0("output/novel_CAGE/CAGE_novel_",NAME,".bed"))
}

annot_writer<-function(LIST,IN){
  NAME<-names(LIST)[IN]
  OUT<-annoGR2DF(LIST[[IN]])
  colnames(OUT)[1]<-"#chr"
  assign(x=NAME,value = OUT)
  write_tsv(x=get(NAME),path = paste0("output/annot_CAGE/CAGE_annot_",NAME,".bed"))
}

for (i in seq_along(CAGE_short_list)){
  CAGE_writer(CAGE_short_list,i,"short")
}
for (i in seq_along(CAGE_long_list)){
  CAGE_writer(CAGE_long_list,i,"long")
}
for (i in seq_along(CAGE_novel_list)){
  novel_writer(CAGE_novel_list,i)
}

for (i in seq_along(CAGE_annot_list)){
  annot_writer(CAGE_annot_list,i)
}

### Intersection of novel beds and WGBS beds were carried out in the terminal using bedtools
# A manual hashed header was added to the WGBS input files as #chr	start	end	name	score	strand
# sort -k1,1 -k2,2n ../input_WGBS/RUMA_CpGs.hypo.dmr.mod > RUMA.bed
# sort -k1,1 -k2,2n ../input_WGBS/MUSLD_CpGs.hypo.dmr.mod > MUSLD.bed
# for b in *.bed;do sed -i 's/ /\t/'g $b;done

#For findign the interecting peaks with hypomethylation peaks 
#bedtools intersect -b ./RUMA.bed -a ../output/novel_CAGE/CAGE_novel_Sk.bed> MUSLD_CAGE_novel_hypo.bed
#bedtools intersect -b ./MUSLD.bed -a ../output/novel_CAGE/CAGE_novel_Skeletal_Muscle_longissimus_dorsi.bed > MUSLD_CAGE_novel_hypo.bed
#bedtools intersect -b ./RUMA.bed -a ../output/annot_CAGE/CAGE_annot_Rumen_atrium.bed > RUMA_CAGE_annot_hypo.bed
#bedtools intersect -b ./MUSLD.bed -a ../output/annot_CAGE/CAGE_annot_Skeletal_Muscle_longissimus_dorsi.bed > MUSLD_CAGE_annot_hypo.bed

#Both novel and annotation regions were validated by WGBS annotated that are hypo and novel that are hypo 
RUMA_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/RUMA_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
RUMA_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/RUMA_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
RUMA_cage_annot_total<-annoGR2DF(CAGE_annot_list$Rumen_atrium)
RUMA_cage_novel_total<-annoGR2DF(CAGE_novel_list$Rumen_atrium)

MUSLD_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/MUSLD_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
MUSLD_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/MUSLD_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
MUSLD_cage_annot_total<-annoGR2DF(CAGE_annot_list$Skeletal_Muscle_longissimus_dorsi)
MUSLD_cage_novel_total<-annoGR2DF(CAGE_novel_list$Skeletal_Muscle_longissimus_dorsi)

AM_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/AM_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
AM_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/AM_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
AM_cage_annot_total<-annoGR2DF(CAGE_annot_list$Alveolar_Macrophages)
AM_cage_novel_total<-annoGR2DF(CAGE_novel_list$Alveolar_Macrophages)

CER_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/CER_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
CER_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/CER_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
CER_cage_annot_total<-annoGR2DF(CAGE_annot_list$Cerebellum)
CER_cage_novel_total<-annoGR2DF(CAGE_novel_list$Cerebellum)

COR_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/COR_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
COR_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/COR_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
COR_cage_annot_total<-annoGR2DF(CAGE_annot_list$Cerebral_Cortex)
COR_cage_novel_total<-annoGR2DF(CAGE_novel_list$Cerebral_Cortex)

LNG_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/LNG_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
LNG_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/LNG_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
LNG_cage_annot_total<-annoGR2DF(CAGE_annot_list$Lung)
LNG_cage_novel_total<-annoGR2DF(CAGE_novel_list$Lung)

OVY_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/OVY_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
OVY_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/OVY_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
OVY_cage_annot_total<-annoGR2DF(CAGE_annot_list$Ovary)
OVY_cage_novel_total<-annoGR2DF(CAGE_novel_list$Ovary)

MUSBF_cage_novel_hypo<-read_tsv("CAGE_WGBS_overlay/MUSBF_CAGE_novel_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
MUSBF_cage_annot_hypo<-read_tsv("CAGE_WGBS_overlay/MUSBF_CAGE_annot_hypo.bed",col_names = c("chr","start","end","width","strand","score"))
MUSBF_cage_annot_total<-annoGR2DF(CAGE_annot_list$Skeletal_Muscle_biceps_femoris)
MUSBF_cage_novel_total<-annoGR2DF(CAGE_novel_list$Skeletal_Muscle_biceps_femoris)


#Verification by WGBS

data.frame(
  Tissue = c("Alveolar_Macrophages",
             "Cerebellum",
             "Cerebral_Cortex",
             "Lung",
             "Ovary",
             "Rumen_atrium",
             "Skeletal_Muscle_biceps_femoris",
             "Skeletal_Muscle_longissimus_dorsi"),
  annot_hypo = c(
    length(AM_cage_annot_hypo$score),
    length(CER_cage_annot_hypo$score),
    length(COR_cage_annot_hypo$score),
    length(LNG_cage_annot_hypo$score),
    length(OVY_cage_annot_hypo$score),
    length(RUMA_cage_annot_hypo$score),
    length(MUSBF_cage_annot_hypo$score),
    length(MUSLD_cage_annot_hypo$score)
  ),
  annot_total = c(
    length(AM_cage_annot_total$score),
    length(CER_cage_annot_total$score),
    length(COR_cage_annot_total$score),
    length(LNG_cage_annot_total$score),
    length(OVY_cage_annot_total$score),
    length(RUMA_cage_annot_total$score),
    length(MUSBF_cage_annot_total$score),
    length(MUSLD_cage_annot_total$score)
  ),
  novel_hypo = c(
    length(AM_cage_novel_hypo$score),
    length(CER_cage_novel_hypo$score),
    length(COR_cage_novel_hypo$score),
    length(LNG_cage_novel_hypo$score),
    length(OVY_cage_novel_hypo$score),
    length(RUMA_cage_novel_hypo$score),
    length(MUSBF_cage_novel_hypo$score),
    length(MUSLD_cage_novel_hypo$score)
  ),
  novel_total = c(
    length(AM_cage_novel_total$score),
    length(CER_cage_novel_total$score),
    length(COR_cage_novel_total$score),
    length(LNG_cage_novel_total$score),
    length(OVY_cage_novel_total$score),
    length(RUMA_cage_novel_total$score),
    length(MUSBF_cage_novel_total$score),
    length(MUSLD_cage_novel_total$score)
  )
) %>%
  mutate(
    annot_novel_hypo = c(annot_hypo + novel_hypo),
    annot_novel_total = c(annot_total + novel_total),
    novel_noise = novel_total - novel_hypo,
    annot_noise = annot_total - annot_hypo
  ) %>%
  reshape2::melt() -> CAGE_WGBS_metrics

  ggplot(CAGE_WGBS_metrics %>% filter(variable!="annot_total",
                                      variable!="novel_total",
                                      variable!="annot_novel_hypo",
                                      variable!="annot_novel_total"),
         aes(Tissue,value,fill=variable,group=variable)) + 
  geom_bar(position="stack",stat="identity") +
  geom_text(aes(label=format(x = value,big.mark = ","),group=variable),position = position_stack(vjust = 0.5)) +
  geom_label(data = CAGE_WGBS_metrics %>% filter(variable =="annot_novel_total"),
             aes(label=paste0("n=",format(x = value,big.mark = ",")),y=30000),show.legend = FALSE)+
  scale_fill_brewer(palette = "RdYlBu",name="Categories",
                    labels=c("Annotated + HypoCpG",
                             "Annotated w/o",
                             "Total (n)",
                             "Novel + HypoCpG",
                             "Novel w/o")) +
    ylab("CAGE tag cluster counts") +
  theme_classic() + ggtitle("CAGE tags overlay with WGBS Hypo-methylation (CpG) sites") -> plot_cage_wgbs

  
plot_cage_wgbs + theme(axis.text.x.bottom = element_text(angle = 30,hjust = 1))

CAGE_WGBS_metrics %>% 
  pivot_wider(id_cols = Tissue,names_from = variable,values_from = value) %>%
  mutate(`Annotated + HypoCpG`=percent(annot_hypo/annot_total,accuracy = 1),
         `Annotated w/o`=percent(annot_noise/annot_total,accuracy = 1),
         `Novel + HypoCpG`=percent(novel_hypo/novel_total,accuracy = 1),
         `Novel w/o`=percent(novel_noise/novel_total,accuracy = 1)) %>%
  select(Tissue,`Annotated + HypoCpG`,`Annotated w/o`,`Novel + HypoCpG`,`Novel w/o`) %>%
  grid.table()

#Matching MI biological distance between tissues in the RNA data

RNA_invert<-t(RNA[-1])
colnames(RNA_invert)<-RNA$txID
RNA_MI<-MIdist(RNA_invert)

cowplot::plot_grid(
factoextra::fviz_dend(
  hclust(RNA_MI, "complete"),
  k = 12,
  k_colors = "lancet",
  type = "phylogenic",
  repel = T,
  color_labels_by_k = T,cex = 0.7,
  phylo_layout = "layout.auto"
),
factoextra::fviz_dend(
  hclust(TSSs_MI, "complete"),
  k = 12,
  k_colors = "lancet",
  type = "phylogenic",
  repel = T,
  color_labels_by_k = T,cex = 0.7,
  phylo_layout = "layout.auto"
),
labels =c("RNA-Seq TPM clusters","Unidirectional TSSs clusters"))

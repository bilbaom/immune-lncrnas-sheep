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


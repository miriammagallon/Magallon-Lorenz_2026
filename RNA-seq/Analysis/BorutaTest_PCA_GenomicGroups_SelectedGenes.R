##################################
# PCA of DEG Genomic.group.order #
##################################
####### Packages Needed
library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(yaml)
library(ggplot2)
library(ggbeeswarm)
library(pheatmap)
library(viridis)
if(!require("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(DOSE)
# library(gage)
library(apeglm)
library(biomaRt)
library(ggrepel)
library(RColorBrewer)
library(gplots)
library(plyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(regioneR)
library(randomForest)
library(caret)
# if(!require("Boruta", quietly = TRUE)) BiocManager::install("Boruta")
library(Boruta)
library(pROC)
#######################  Loading Functions #########################
source(file= "./Analysis/rna_seq_Functions.R")

#PCAplot deseq function modified
plotPCA <- function (object, intgroup = "condition", pcx1 = "PC1", pcx2="PC2", ntop = 500, returnData = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  pcax <- pca$x
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                               drop = FALSE])
  if (length(intgroup) > 1) {
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }else {
    group <- colData(object)[[intgroup]]
  }
  
  d <- cbind(pcax, group,
             intgroup.df, name = colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:ncol(pcax)]
    return(d)
  }
  ggplot(data = d, aes_string(x = pcx1, y = pcx2, color = "group")) +
    geom_point(size = 3) + xlab(paste0(pcx1,": ", round(percentVar[1] *
                                                          100), "% variance")) + ylab(paste0(pcx2,": ", round(percentVar[2] *
                                                                                                                100), "% variance")) + coord_fixed()
}

#color ramp
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#######################     Parameters     ##########################

#### loading the sample data information
sample.data <- read.table(file = "./Sample.info.csv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

#Our set and Miller set
sample.data <- sample.data[sample.data$Real.ValSet %in% c("Real","ValSet"),]
sample.data <- sample.data[sample.data$Analysis!="",]
sample.data <- sample.data[!grepl("ChiPing",sample.data$Sample.name),] 

sample.group <- "Genomic.group.order"
sample.data$Genomic.group.order <- as.character(sample.data$Genomic.group.order)
sample.data$Genomic.group.order[is.na(sample.data$Genomic.group.order )] <- ""


if(sample.group == "Genomic.group.order"){
  sample.type <- unique(sample.data$Genomic.group.order)
  sample.type <- sample.type[sample.type!=""]
}


###### Parameters
sample.data$Sample.name[sample.data$Sample.name =="MPNST_SP_09"] <- "15PNF"
file.names <- sample.data$Sample.name
sample.names <- sample.data$Sample.name
sample.type <-unique(sample.data[,sample.group][sample.data[,sample.group]!=""])

#Salmon parameters
output.suffix <- "_quant"
output.quants <- "./Salmon"
# threads <- 8

# Tximport paramenters
orgdb <- org.Hs.eg.db
org.columns <- "SYMBOL"
org.keytype <- "REFSEQ"

# DESeq2 parameters
filt.min.reads <- 5
filt.min.samples <- 3 #chaged the number of samples as minimun we had 5 replic

#Enrichment parameters
goterms <- c("BP", "CC", "MF","KEGG")

######### Directories
res.dir <- file.path("results",sample.group)
if(!file.exists(res.dir)) dir.create(res.dir)

pca.dir <-  file.path(res.dir,  "PCAplot")
if(!file.exists(pca.dir)) dir.create(pca.dir)

distance.dir <- file.path(res.dir,  "DistancePlots")
if(!file.exists(distance.dir)) dir.create(distance.dir)

heatmap.QC.dir <- file.path(res.dir, "heatmap.QC")
if(!file.exists(heatmap.QC.dir)) dir.create(heatmap.QC.dir)

heatmap.dir <- file.path(res.dir,"DEG_heatmap")
if(!file.exists(heatmap.dir)) dir.create(heatmap.dir)

DEG.dir <- file.path(res.dir, "DESeqResults")
if(!file.exists(DEG.dir)) dir.create(DEG.dir)

selected.genes.dir <- file.path(res.dir, "DEGgenes")
if(!file.exists(selected.genes.dir)) dir.create(selected.genes.dir)

markers.dir <- file.path(res.dir, "markers")
if(!file.exists(markers.dir)) dir.create(markers.dir)

gsea.dir <- file.path(res.dir,  "GSEAGO")
if(!file.exists(gsea.dir)) dir.create(gsea.dir)
for(term in goterms){
  if(!file.exists(file.path(gsea.dir, term))) dir.create(file.path(gsea.dir, term))
}

annot.dir <- file.path(res.dir, "Annotation")
if(!file.exists(annot.dir)) dir.create(annot.dir)

gene.expression.dir <- file.path(res.dir,  "GeneExpressionPlot")
if(!file.exists(gene.expression.dir)) dir.create(gene.expression.dir)

cluster.dir <- file.path(res.dir,"GeneCluster")
if(!file.exists(cluster.dir)) dir.create(cluster.dir)


################## Color of samples ################
annot.column <- data.frame(sample.name = sample.data$Sample.type,
                           SpvsNF1 =sample.data$SpvsNF1,
                           TumorvsCL = sample.data$TumorvsCL,
                           BeningvsMalign = sample.data$BenignvsMalign,
                           Purity = sample.data$Tumor.purity,
                           TSG.inact = sample.data$TSG_inactivation,
                           NF1.inact = sample.data$NF1,
                           CDKN2A.inact = sample.data$CDKN2A,
                           CDKN2AH.inact = sample.data$CDKN2A.SV.hotspot,
                           PRC2.inact = sample.data$PRC2,
                           H3K27me3 = sample.data$H3K27me3,
                           TP53.inact = sample.data$TP53,
                           PTEN.inact = sample.data$PTEN,
                           RB1.inact = sample.data$RB1,
                           CN.LOH.regions = sample.data$CN.LOH.0.5,
                           LOH.prop = sample.data$LOH.0.5,
                           Activating.mut = sample.data$Activating.mut,
                           Mut.freq = sample.data$Mut.freq,
                           Even.chr = sample.data$Even.chr.50,
                           platform = sample.data$Seq.Platform,
                           Methylation.id =sample.data$Methylome_ID,
                           Genomic.group = sample.data$Genomic.group.order,
                           Cluster = sample.data$Clusters,
                           chr1.unbalance = sample.data$chr1.unbalance)

rownames(annot.column) <- sample.data$Sample.name
#Samplet type colors
sample.name.col <- unique(sample.data$Group.color)
names(sample.name.col) <- unique(sample.data$Sample.type)
sample.name.col <- unique(sample.name.col)
names(sample.name.col) <- unique(sample.data$Sample.type)

#Platform colors
col.paltform <- c(ILL = "yellow", BGI ="orange")

#Purity color ramp
purity <- annot.column$Purity
pal <- colorRampPalette(c("darkblue","white"))
color.purity <- map2color(x = purity,pal =pal(15) ,limits = c(0,1))
names(color.purity) <-unique(sample.data$Sample.name)
color.purity[is.na(color.purity)] <- "grey"

#Mut.freq
mut.freq <- annot.column$Mut.freq
pal <- colorRampPalette(c("white","orange"))
color.mut.freq <- map2color(x = mut.freq,pal =pal(19) ,limits = c(0,3e5))
names(color.mut.freq) <-unique(sample.data$Sample.name)
color.mut.freq[is.na(color.mut.freq)] <- "grey"

# TSGs color
color.tsg <- as.character(annot.column$TSG.inact)
color.tsg[is.na(color.tsg)] <- "grey93"
color.tsg[color.tsg == "FALSE"] <- "grey72"
color.tsg[color.tsg == "TRUE"] <- "#FF6625"
color.tsg <- unique(color.tsg)
if(length(color.tsg)==3){
  names(color.tsg) <- c("FALSE","TRUE","")
}else{
  names(color.tsg) <- c("FALSE","TRUE")
  
}

#TP53 inact color
# color.tp53 <- c("FALSE" = "#73D216FF", "TRUE"="#75507BFF")
color.tp53 <- as.character(annot.column$TP53.inact)
color.tp53[is.na(color.tp53)] <- "grey93"
color.tp53[color.tp53 == "FALSE"] <- "grey72"
color.tp53[color.tp53 == "TRUE"] <- "#D8BDF9"
color.tp53 <- unique(color.tp53)

if(length(color.tp53)==3){
  names(color.tp53) <- c("FALSE","TRUE","")
}else{
  names(color.tp53) <- c("FALSE","TRUE")
  
}

#Methylation group color
color.methyl <- c("white","#572B54","black","#FB6E52","yellow","pink","lightblue","grey","darkblue","darkred","brown")
names(color.methyl) <- unique(sample.data$Methylome_ID)

#NF1/SP color group
color.nf1sp <- c("NF1"= "#5AC1D1", "SP"="#815AD1", "white")

#PRC2inact
color.prc <- c( "TRUE"="#D8BDF9","FALSE" = "grey72","grey93")
# annot.column$PRC2.inact
#Tissuelvs cell line
color.tissuevscl <- c("Tumor"="#D15A80", CL="#B9EDDC")

#Benign vs Malignant
color.bm <- c( Benign = "#4C9B6B", Malignant = "#8A2450")

#NF1 inact
color.nf1 <- as.character(annot.column$NF1.inact)
color.nf1[is.na(color.nf1)] <- "grey93"
color.nf1[color.nf1 == "FALSE"] <- "grey72"
color.nf1[color.nf1 == "TRUE"] <- "#D8BDF9"
color.nf1 <- unique(color.nf1)
if(length(color.nf1)==3){
  names(color.nf1) <- c("TRUE","FALSE","")
}else{
  names(color.nf1) <- c("TRUE","FALSE")
  
}


#CDKN2A Hotspot
color.cdkn2aH <- as.character(annot.column$CDKN2AH.inact)
color.cdkn2aH[is.na(color.cdkn2aH)] <- "grey93"
color.cdkn2aH[color.cdkn2aH == "FALSE"] <- "grey72"
color.cdkn2aH[color.cdkn2aH == "TRUE"] <- "#8E3DEE"
color.cdkn2aH <- unique(color.cdkn2aH)
# names(color.cdkn2aH) <- c("FALSE","","TRUE")
if(length(color.cdkn2aH)==3){
  names(color.cdkn2aH) <- c("FALSE","","TRUE")
}else{
  names(color.cdkn2aH) <- c("FALSE","TRUE")
  
}

color.cdkn2a <- as.character(annot.column$CDKN2A.inact)
color.cdkn2a[is.na(color.cdkn2a)] <- "grey93"
color.cdkn2a[color.cdkn2a == "FALSE"] <- "grey72"
color.cdkn2a[color.cdkn2a == "TRUE"] <- "#D8BDF9"
color.cdkn2a <- unique(color.cdkn2a)
if(length(color.cdkn2a)==3){
  names(color.cdkn2a) <- c("FALSE","","TRUE")
}else{
  names(color.cdkn2a) <- c("FALSE","TRUE")
  
}
#Activating mut
color.actmut <- as.character(annot.column$Activating.mut)
color.actmut[is.na(color.actmut)] <- "grey93"
color.actmut[color.actmut == "FALSE"] <- "grey72"
color.actmut[color.actmut == "TRUE"] <- "#0B007A"
color.actmut <- unique(color.actmut)
# names(color.actmut) <- c("FALSE","TRUE","")
if(length(color.actmut)==3){
  names(color.actmut) <- c("FALSE","TRUE","")
}else{
  names(color.actmut) <- c("FALSE","TRUE")
}

# Genomic groups
# col.groups <- c("#DA6220",  "tomato","goldenrod","#C55454FF","#E09B9B","white")
# names(col.groups) <-c("1","2","3","4","5","")
col.groups <- c("#7C2C00",  "tomato","goldenrod","#E09B9B","grey93")
names(col.groups) <-c("1","2","3","4","")
# H3K27me3
col.H3k27me3 <- as.character(annot.column$H3K27me3)
col.H3k27me3[col.H3k27me3 ==""] <- "grey93"
col.H3k27me3[col.H3k27me3 %in% "Loss"] <- "#52B5FB"
col.H3k27me3[col.H3k27me3  %in% "Retained"] <- "#FB5272"
col.H3k27me3 <- unique(col.H3k27me3)
if(length(col.H3k27me3)>1){
  names(col.H3k27me3)<- c("","Loss","Retained")
} else{
  names(col.H3k27me3) <- c("")
  
}

# >50% CN-LOH regions
color.cn.loh <- c( "TRUE"="#CD853F","FALSE" = "grey72","grey93")

# >50%  Whole genome LOH
color.wgs.loh <- c( "TRUE"="#6591A5","FALSE" = "grey72","grey93")

# even chromosomes
color.even.chr <- c( "TRUE"="#73D216FF","FALSE" = "grey72", "grey93")

# chr1 unbalance
color.chr1.unb <- c( "TRUE"="#87CEFA","FALSE" = "grey72", "grey93")

#clusters
col.clusters <- c(C1= "#7BA7FFFF", C2 = "#FFEC6DFF","Grey"= "grey","grey92")

annot.column$LOH.prop[annot.column$LOH.prop == ""] <- "FALSE"
annot.column$Even.chr[annot.column$Even.chr == ""] <- "TRUE"
mycolors <- list(sample.name = sample.name.col,
                 SpvsNF1 = color.nf1sp,
                 platform = col.paltform,
                 TumorvsCL= color.tissuevscl,
                 Purity = color.purity,
                 BeningvsMalign = color.bm,
                 NF1.inact =color.nf1,
                 CDKN2A.inact =color.cdkn2a,
                 CDKN2AH.inact =color.cdkn2aH,
                 PRC2.inact = color.prc,
                 TSG.inact = color.tsg,
                 TP53.inact = color.tp53,
                 PTEN.inact = color.tp53,
                 RB1.inact = color.tp53,
                 Activating.mut = color.actmut,
                 CN.LOH.regions = color.cn.loh,
                 LOH.prop = color.wgs.loh,
                 Mut.freq = color.mut.freq,
                 Even.chr = color.even.chr,
                 Methylation.id = color.methyl,
                 Genomic.group = col.groups,
                 H3K27me3 = col.H3k27me3,
                 Cluster = col.clusters,
                 chr1.unbalance = color.chr1.unb
                 
)
pvalue <- 0.05

annot.column[sapply(annot.column, is.logical)] <- lapply(annot.column[sapply(annot.column,is.logical)],as.character)
annot.column[sapply(annot.column, is.character)][sapply(annot.column[sapply(annot.column, is.character)],is.na)] <-""


##################              Tximport              #####################
salmonquants.fl <-file.path(output.quants, paste0(file.names,"_quant/","quant.sf"))
names(salmonquants.fl)<- file.names
salmonquants.fl <- salmonquants.fl[file.exists(salmonquants.fl)] 
txi.salmon <- importQuantsData(quant.files = salmonquants.fl, orgdb = orgdb,
                               orgdb.keytype = orgdb.keytype, org.columns = org.columns,
                               tximpot.type = "salmon")
#Change Sample name
colnames(txi.salmon$abundance) <- names(salmonquants.fl)
colnames(txi.salmon$counts) <- names(salmonquants.fl)
colnames(txi.salmon$length) <- names(salmonquants.fl)
# txi.salmon <- selectDataFromTximport(tximport = txi.salmon, sample_names = sample.names[-c(12,15)])

##################          Quality Control       ######################
sample.data <-sample.data[sample.data$Sample.name %in% names(salmonquants.fl),]
# rownames(sample.data) <- sample.data$Sample.name
# We import the data from tximport output to analyse it using DESeq2
dds.txi <- DESeqDataSetFromTximport(txi = txi.salmon, colData = sample.data, formula(paste0("~", sample.group)))
dim(dds.txi)

# Estimation of Size Factors
dds <- estimateSizeFactors(dds.txi)

# Filtering rawdata: removing data without expression
keep <- rowSums(counts(dds)>=filt.min.reads)>=filt.min.samples
dds <- dds[keep,]

#Normalization
dds.deg <- DESeq(object = dds, test = "Wald")
rdata.file.rlog <- file.path(res.dir,"dds.rlog_AllSamplesAllGenes_RealSetValSet.RData")
if(file.exists(rdata.file.rlog)){
  load(rdata.file.rlog)
}
#################      Differential gene expression      ########################
#Sample.data 
sample.sel.data <- sample.data[sample.data[,sample.group]!="", ]
sample.sel.data <- sample.sel.data[grepl("Miller",sample.sel.data$Sample.name),]
sample.sel.data <- sample.sel.data[sample.sel.data$Genomic.group.order %in% c(1,2,3,4),]

sample.sel.data$Genomic.group.order <- factor(sample.sel.data$Genomic.group.order)

###### Tximport ########
#File names
file.names <- sample.sel.data$Sample.name

salmonquants.fl <-file.path(output.quants, paste0(file.names,"_quant/","quant.sf"))
names(salmonquants.fl)<- file.names
salmonquants.fl <- salmonquants.fl[file.exists(salmonquants.fl)] 
txi.sel<- importQuantsData(quant.files = salmonquants.fl, orgdb = orgdb,
                           orgdb.keytype = orgdb.keytype, org.columns = org.columns,
                           tximpot.type = "salmon")
#Change Sample name
colnames(txi.sel$abundance) <- names(salmonquants.fl)
colnames(txi.sel$counts) <- names(salmonquants.fl)
colnames(txi.sel$length) <- names(salmonquants.fl)

#Filtering deseqData
filtered.dds <- getFilteredDDS(tximport = txi.sel, samples_group = sample.group, samples_df = sample.sel.data, filter_min_reads = filt.min.reads, filter_min_samples = filt.min.samples)

dds.rlogsel <-rlog(filtered.dds)

################ PCA representation ###################

###### Load DEG lists after groups comparison #######

#G1G2vsG3G4 up in G1G2
genes.results <- read.table(file.path(res.dir,"DEGgenes","ResultTable_G1G2vsG3G4_Genomic.group.order_G1G2vsG3G4_ValSet.csv"), sep = "\t",header = T)
genes.G1G2 <- genes.results$genes[genes.results$log2FoldChange<(0)]
genes.G1vsG2 <- read.table(file.path(res.dir,"DEGgenes","ResultTable_1vs2_Genomic.group.order_ValSet.csv"), sep = "\t",header = T)
genes.G1vsG2 <- genes.G1vsG2$genes


genes.G1vsG3 <- read.table(file.path(res.dir,"DEGgenes","ResultTable_1vs3_Genomic.group.order_ValSet.csv"), sep = "\t",header = T)
genes.G2vsG3 <- read.table(file.path(res.dir,"DEGgenes","ResultTable_2vs3_Genomic.group.order_ValSet.csv"), sep = "\t",header = T)
genes.G3vsG4 <- read.table(file.path(res.dir,"DEGgenes","ResultTable_4vs3_Genomic.group.order_ValSet.csv"), sep = "\t",header = T)
genes.G3 <- genes.G1vsG3$genes[genes.G1vsG3$genes[genes.G1vsG3$log2FoldChange>0]%in%genes.G2vsG3$genes[genes.G2vsG3$log2FoldChange>0]]
genes.G3 <- genes.G3[genes.G3%in% genes.G3vsG4$genes[genes.G3vsG4$log2FoldChange>0]]

genes.G4 <- genes.G3vsG4$genes[genes.G3vsG4$log2FoldChange<(0)]


#from this lists we perform a random forest for selecting the most influencing genes
genes <- unique(c(genes.G1G2, genes.G1vsG2, genes.G3, genes.G4)) 
length(genes)

read.table(file.path("./DataFigures/1703genes_forBorutaTest.csv"),sep = "\t",header = FALSE )



#INTENTANDO COGER OTRO NUMERO DE GENEES
# genes <- unique(c(genes.G1G2, genes.results$genes, genes.results.1$genes, genes.results.2$genes, genes.results.3$genes,
#                   genes.results.4$genes,
#                   genes.results.5$genes,
#                   genes.results.6$genes))
# length(genes)
# "SOX9"%in%genes

# genes <- geness[geness%in%genes]
# table(genes %in% "CDKN2A")
# genes <- c(genes,"CDKN2A")
# genes <- unique(c("CLEC2A","VGLL2","RBM24","SCN3B","NEB","XAGE1A","CDKN2A","NFGR","GRIA4","INA","GRIK1",
#            "ZIC3","CA13","NXPH2","SP8","GPC5","HOXD13","ETV4","PPEF1",genes.g3))

# enrich <- getAllEnrichTermList(genes.G3)
# enrich$BP$geneID

############ Gene selection using Boruta's statistic and  Recursive  Feature Elimination with Cross-Validation (RFECV)  ###########

# This two-step approach is powerful:
# 1- Step 1: Boruta: we use Boruta on all 1795 genes to filter out the "noise and identify the entire set of "all-relevant" genes.
#            This gives a robust, reproducible, and biologically meaningfull list
#
# 2- Step 2: RFECV: we take the smaller list of "Confirmed" genes from Boruta and use RFECV with CV on it.
#            This finds the "minimal-optimal" subset for the best predictive performance


length(genes)
gene.exp <- assay(dds.rlogsel)[rownames(dds.rlogsel)%in% genes,] # Only Miller Samples

# gene.exp <- assay(dds.rlog)[rownames(dds.rlog)%in%genes, colData(dds.rlog)$Genomic.group.order %in% c(1,2,3,4)]

nrow(gene.exp)
sample.groups <- colData(dds.rlogsel)[,sample.group]
# sample.groups <- colData(dds.rlog)[colData(dds.rlog)$Genomic.group.order %in% c(1,2,3,4),sample.group]
# sample.groups <- factor(as.character(sample.groups), levels =c(1,2,3,4))
expr_df <- as.data.frame(t(gene.exp))
colnames(expr_df) <- make.names(colnames(expr_df)) # fix names for formula
# expr_df$Group  <- sample.groups

set.seed(123)
boruta.r <- Boruta(x = expr_df, y = sample.groups, doTrace =2, maxRuns = 5000,seeds = 123)
boruta.r <- Boruta::TentativeRoughFix(boruta.r)
boruta.genes <- getSelectedAttributes(boruta.r, withTentative =F)
# n.runs <- 100
# all.boruta.results <- vector(mode ="list", length = n.runs)
# for(i in 1:n.runs){
#   set.seed(i)
#   set.seed(123)
#   boruta.r <- Boruta(x = expr_df, y = sample.groups, doTrace =0, maxRuns = 5000)
#   boruta.genes <- getSelectedAttributes(boruta.r, withTentative =F)
#   all.boruta.results[[i]] <-boruta.genes
#   if(i%%10 ==0){print(paste("Completed run", i,"of", n.runs))}
# }
# save(all.boruta.results, file = file.path(markers.dir,"Boruta_100times_genes_NoTentative.RData") )
# load(file.path(markers.dir,"Boruta_100times_genes_NoTentative.RData"))
# all.boruta.results
# #ulist all 100 list into one giant vector
# all.genes.vector <- unlist(all.boruta.results)
# # Count the frequency of each gene
# gene.freq <- table(all.genes.vector)
# 
# # Sort by frequency
# stable.genes.tb <- sort(gene.freq, decreasing = T)
# 
# # we stablish a thr to select genes that appeared in at least 95% of runs
# thr <- 95
# 
# genes.sel <- names(stable.genes.tb[stable.genes.tb>=thr])
confirmed.genes.b <- expr_df[,boruta.genes]

# Get your reproducible, Biologically Relevan List
# boruta.genes <- getSelectedAttributes(boruta.r, withTentative =T)
# plot(boruta.r, xlab = "Genes", ylab = "Importance")
message("These are your biologically relevant genes")
# print(genes.sel)
genes.sel <- gsub(pattern = "\\.",replacement = "-",boruta.genes)

# Define the "control settings for RFECV
rfe_control <- rfeControl(functions = rfFuncs, method = "repeatedcv", number = 10, repeats = 5)

# Run the RFECV Algorithm using the Boruta-reduced data: A simple backwards selection--> Recursive feature elimination algorithm
set.seed(123)
rfe_result <- rfe(x = confirmed.genes.b, y = sample.groups, 
                  sizes=c(1:10,15:20,40,60), # Which feature set sizes to test
                  rfeControl = rfe_control)

#Get the final, smallest, optimal list of genes
final_signature <- predictors (rfe_result)
genes.sel <- gsub(pattern = "\\.",replacement = "-",final_signature)

# This plot shows you the model accuracy at each gene set size
plot(rfe_result, type =c("g","o"))
# genes.sel <- head(genes.sel,30)
length(genes.sel)


# write.table(genes.sel,file.path(markers.dir,"Boruta_RFECV_Confirmed_FinalSignature_OrderMostContibutingGenes_76genes_5000Runs.csv"),col.names = F,row.names = F,sep = "\t",quote = F)
# write.table(genes.sel,file.path(markers.dir,"Boruta_RFECV_Confirmed_FinalSignature_OrderMostContibutingGenes_76genes_5000Runs.csv"),col.names = F,row.names = F,sep = "\t",quote = F)

genes.sel <- read.table(file.path(markers.dir,"Boruta_RFECV_Confirmed_FinalSignature_OrderMostContibutingGenes_76genes_5000Runs.csv"), header = F, sep = "\t")[,1]

###### PCA representation #######
genes.sel <- c("CCNE1","PRR4", "DMRTA2","MYOG","LGR5", "BARX2","KLRC2","EMX1","INA")
genes.sel <- c("MYOD1","PAX7","DES","MYOG","PAX3","FAP")

# genes.sel <- ("APOD","DCN")
# genes.sel <- c("PRR4","CLEC2A","CCNE1","KLRC2","FEZF1","ZIC3","SALL3","NKX2-4")
dds.rlogs <- dds.rlog[rownames(dds.rlog)%in%genes.sel, colData(dds.rlog)$Genomic.group.order %in% c(1,2,3,4)]
# dds.rlogs <- dds.rlogsel[rownames(dds.rlogsel)%in%genes.sel, colData(dds.rlogsel)$Genomic.group.order %in% c(1,2,3,4)]
# dds.rlogs <- dds.rlogs[,colnames(dds.rlogs)!="MPNST-32_Miller"]
pca.pseudocounts <- dds.rlogs

#Intgroups are the colnames present in colData(dds.rlogs)
intgroup <- c("Genomic.group.order")
# intgroup <- c("LOH.0.5")
# intgroup <- c("RB1")

# colData(pca.pseudocounts)
pcaData <- plotPCA(pca.pseudocounts, intgroup=intgroup, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
color.sample <- mycolors[["Genomic.group"]]
# color.sample <- mycolors[["LOH.prop"]]


pp <- ggplot2::ggplot(pcaData, aes(PC1, PC2, color = pcaData[,intgroup] ))+
  geom_point(size = 5, stroke =3,pch = 16)+ 
  scale_color_manual(name = intgroup, values =color.sample)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  # ggrepel::geom_label_repel(aes(label = as.character(rownames(pcaData))), show.legend = FALSE, size = 10)+
  theme(legend.title = element_text(size=60), 
        legend.text = element_text( size=60), 
        legend.position = "top",
        axis.text=element_text(size=70),
        axis.title=element_text(size=70),
        plot.title = element_text(size=50),
        panel.background = element_rect(fill="white"), 
        panel.grid.major = element_line(colour = "grey50"), plot.margin =ggplot2::margin(t=30,r=100,b=30,l=30))+
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 15)),
         shape = guide_legend( override.aes = list(size = 5, color ="#7C8e8F", fill = "#7C8e8F" )))

# ggsave(plot = pp,filename =file.path(pca.dir, paste0("PCA_PC1-PC2_TFgenes_AllSamples.", intgroup,".", sample.group,".nonames_IntBiol_RealSet_v4.svg")), width =40, height = 28)
ggsave(plot = pp, filename =file.path(pca.dir, paste0("PCA_PC1-PC2_SelectedGenes_Boruta_RFECV_Confirmed_FinalSignature_OrderMostContibutingGenes_84genes_5000Runs_87Accuracy.", intgroup,".", sample.group,".names_IntBiol_RealValset.svg")), width =40, height = 28)

###### Heatmap representation #####
genes.sel <- c("PTEN","CDKN2A","NGFR")
#Selection of genes to represent
# data_subset_norm <- dds.rlogsel[,colData(dds.rlogsel)$Genomic.group.order %in%c(1:4)] #Miller
data_subset_norm <- dds.rlog[,colData(dds.rlog)$Genomic.group.order %in%c(1:4)] #All Tumors
# data_subset_norm <- dds.rlog[,grepl("MPNST",colData(dds.rlog)$Sample.type)] #All Tumors and Cell lines

#Change Sample name for the paper
colnames(colData(data_subset_norm))
paper.names <- colData(data_subset_norm)$Paper.names
# colnames(data_subset_norm) <- paper.names
# data_subset_norm <- dds.rlogsel[,colData(dds.rlogsel)$Genomic.group.order %in%c("G1G2","G3G4")]

data_subset_norm <- assay(data_subset_norm)[rownames(data_subset_norm) %in%genes.sel,]
data_subset_norm <- t(apply(data_subset_norm, 1, cal_z_score))

#Column order
anc <- annot.column[rownames(annot.column) %in% colnames(data_subset_norm),] 
rownames(anc) <- paper.names
colnames(data_subset_norm) <- paper.names
anc$Genomic.group.order <- factor(anc$Genomic.group,levels =c(1,2,3,4,""))
anc <- anc[order(anc$Genomic.group),]
sn.order <- rownames(anc)
names(sn.order)
data_subset_norm <- data_subset_norm[,sn.order] 

# Column annotation 
# colnames(anc)
# anno.col <- anc[,c(1,6,8,10,11,12,2,14,20)] 
anno.col <- anc[,c(11,14,12,10,8,6,2,22)] 
anno.col <- anc[,c(20,11)] 

# png(filename = file.path(heatmap.dir,paste0("DEG_G1vsG2_",sample.group,"_ColumnUNSupervised_IntBiol_RealValset.v2.png")), width = 1200, height = 1000)
png(filename = file.path(heatmap.dir,paste0("Boruta_RFECV_Confirmed_FinalSignature_OrderMostContibutingGenes_84genes_5000Runs_87Accuracy_",sample.group,"_ColumnUNSupervised_IntBiol_RealValset.v2.png")), width = 1200, height = 1000)

# png(filename = file.path(heatmap.dir,paste0("30PC3genes_supervised.2_AllSamples.IntBiol.png")), width = 1500, height = 2500)
# png(filename = file.path(heatmap.dir,paste0("MesenchymalMarkers_Itzi_",sample.group,"_UnSupervised.AllSamples.Realset.png")), width = 1200, height = 1000)

pheatmap(data_subset_norm,
         # annotation_row = anno.row,
         annotation_col = anno.col,
         annotation_colors = mycolors,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows = T,
         # breaks = myBreaks,
         fontsize = 10, fontsize_row = 10,
         color = bluered(80),
         show_rownames =T,
         annotation_legend = TRUE, 
         legend = TRUE,
         annotation_names_row = F,
         annotation_names_col = T,
         treeheight_row = 50, 
         treeheight_col = 50,
         scale ="row",
         lwd = 1)
dev.off() 


####### Gene expression of selected Genes ######
##### Boxplot and barplot representation ####
genes.sel <- c("PYCARD","IL1B","NLRP3","IL18","ZIC3")
genes.sel <- c("SNAI1","ZEB1","NRG1","TGFB2")
genes.sel <- c("PTEN","RB1","SOX9","SOX10","CDH19","S100B","PLP1","TWIST1","PRRX1","NGFR")
genes.sel <- c("NTRK3","NTRK2","MYOD1","NTRK1")
genomic.groups.data <-sample.data[sample.data$Genomic.group.order %in% c(1,2,3,4),]
#Only Miller samples
genomic.groups.data <- genomic.groups.data[grepl("Miller",genomic.groups.data$Sample.name),]
txi.sel <- selectDataFromTximport(tximport = txi.salmon,sample_names = genomic.groups.data$Sample.name )
filtered.dds.groups <- getFilteredDDS(tximport = txi.sel, samples_group = sample.group, samples_df = genomic.groups.data, filter_min_reads = filt.min.reads, filter_min_samples = filt.min.samples)

filt.counts <- counts(filtered.dds.groups, normalize = TRUE)
write.table(filt.counts,file = file.path("/imppc/labs/eslab/mmagallon/Projects/RNA-seq_MesenchimalDiff/results/Genomic.group.order/rna_se_expression_MillerSamples.csv"),sep="\t",col.names = T, row.names = T)
#Selecting genes
filt.counts <- filt.counts[rownames(filt.counts) %in% c(genes.sel),]
fct <- t(filt.counts)
fct <- data.frame(fct)

#Adding paper names
rownames(fct) <- colData(filtered.dds.groups)$Paper.names
colnames(filt.counts) <-rownames(fct)
stages <- factor(genomic.groups.data$Genomic.group.order, levels = c(1,2,3,4))
fct<- cbind(fct, stages)
fct <- fct[order(fct$stages),]

g.group.col <- mycolors$Genomic.group[-length(mycolors$Genomic.group)]
color.data <- as.character(fct$stages)
for(i in seq_len(length(g.group.col))){
  col.g <- g.group.col[i]
  color.data[color.data %in% names(col.g)] <- col.g
  # names(color.data[color.data %in% names(col.g)]) <- names(col.g)
}
i=2
#Boxplot representation

library(plyr) # Assuming you are using plyr for round_any

for(i in seq_len(ncol(fct)-1)){
  
  # 1. DEFINE YOUR ZOOM LIMIT
  # You must change this to your desired zoomed cutoff. 
  # If you use max(fct[,i]), no outliers will be cut off.
  # Example: y_zoom_limit <- 1000  (or maybe a 95th percentile calculation)
  if(names(fct)[i] == "NGFR"){
    y_zoom_limit <- 65000
  } else if(names(fct)[i] == "S100B"){
    y_zoom_limit <- 10000 # <--- UPDATE THIS VALUE
  }else if(names(fct)[i] == "CDH19"){
    y_zoom_limit <- 22000 # <--- UPDATE THIS VALUE
  }else{
    y_zoom_limit <-max(fct[,i])
  }

  
  # Update your 'm' logic to use the zoomed limit instead of the absolute max
  m <- y_zoom_limit
  
  if(m<=999){
    plus <- 50
  } else{
    plus <- 500
  }
  # The axis boundary
  top_limit <- m + plus
  
  pos <- round_any(seq(0, top_limit, by = plus), plus, f = ceiling)
  pos.at  <- c(min(pos), mean(pos), max(pos))
  
  svg(filename = file.path(markers.dir,
                           paste0("Boxplot_", colnames(fct)[i], "_NOnormalizedCounts_GenomicGroupSamples_RF_borut_selection_RealValset_paperNames_76GeneSel_V2.svg")), width = 12, height = 10)
  # png(filename = file.path(markers.dir,
  #                          paste0("Boxplot_", colnames(fct)[i], "_NOnormalizedCounts_GenomicGroupSamples_RF_borut_selection_RealValset_paperNames_1703_76GeneSel_V2.png")), width = 1200, height = 1000)
  # 
  par(mar = c(7,15,7,20)+0.1, mgp = c(10,1,0), yaxs = "r")
  
  # 2. DRAW BOXPLOT WITH ZOOMED YLIM
  bp <- boxplot(fct[,i]~fct$stage, 
                lwd=3,
                xlab = "",
                ylab="",
                main = colnames(fct)[i],
                boxfill=unique(color.data),
                frame=FALSE,
                cex.axis = 5,
                cex.main =3,
                cex.lab = 20,
                xaxt="n", yaxt="n", 
                ylim= c(0, top_limit)) # Uses zoomed limit
  
  # --- 3. UPDATED OUTLIER SQUISH LOGIC + LABELS ---
  # Check if any outliers exceed the boundary
  out_to_squish <- bp$out > top_limit
  
  if(any(out_to_squish)) {
    
    # 3a. Extract the data for only the out-of-bounds outliers
    outlier_values <- bp$out[out_to_squish]
    outlier_groups <- bp$group[out_to_squish]
    
    # 3b. Draw the markers (triangles work better at the plot boundary)
    points(x = outlier_groups,
           y = rep(top_limit, length(outlier_values)), # Plot exactly on the boundary
           pch = 16, col = "red", 
           cex = 1.5) # Match massive axis text sizes
    
    # 3c. ADD THE VALUES (THE MODIFICATION)
    # Using round(..., 0) assuming these are counts.
    text(x = outlier_groups, 
         y = rep(top_limit, length(outlier_values)), # Target coordinate (top of plot)
         labels = paste0(" ", round(outlier_values, 0)), # Rounding for neatness
         col = "red", 
         cex = 1.5, # Large size, matching your publication style (Axis is 5, Legend is 3)
         pos = 1,   # Placed BELOW the coordinate so it's inside the plot area
         offset = 0.8) # Adjust vertical distance from the marker
  }
  
  axis(1, at = c(0:15),labels = rep("",16),lwd.ticks = 0, lwd = 3)
  axis(2, at = pos.at, las=2, lwd = 3, cex.axis=5, ylim =c(0, m + plus))
  
  legend("topright", inset=c(-0.4,0), box.lwd = "", xpd = TRUE, fill = unique(color.data), legend = levels(fct$stage), horiz = F, cex = 3, col = 2)
  
  dev.off()
}




png(filename = file.path(markers.dir,
                         paste0("Boxplot_", colnames(fct)[i], "_NOnormalizedCounts_GenomicGroupSamples_RF_borut_selection_Valset_paperNames_1703_76GeneSel_V2.png")), width = 1200, height = 1000)


###### Barplots

filt.counts <- filt.counts[,rownames(fct)]
i=8
for( i in seq_len(nrow(filt.counts))){
  png(filename = file.path(file.path(markers.dir, paste0("Barplot_",rownames(filt.counts)[i],"_GenomicGroupSamples_Boruta_RFCV_selectedgenes_RealValset_paperNames_1703_76GeneSel_V2_PaperNames.png"))),
      width = 1500, height = 1000)
  par( mar =c(12,10,7,5)+0.1, mgp =c(10,1,0), yaxs = "r")
  barplot(filt.counts[i,], las=2, col = color.data, main= rownames(filt.counts)[i], cex.main = 3, cex.axis = 3, cex.names = 1.5)
  dev.off()
}





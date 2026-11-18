#####################################################
#    Integrative biology: PCAs OurData              #
#####################################################
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

#real set
sample.data <- sample.data[sample.data$Real.ValSet == "Real",]
sample.data <- sample.data[sample.data$Analysis!="",]

sample.group <- "Clusters"
sample.data$Genomic.group.order <- as.character(sample.data$Genomic.group.order)
sample.data$Genomic.group.order[is.na(sample.data$Genomic.group.order )] <- ""


if(sample.group == "Clusters"){
  sample.type <- unique(sample.data$Clusters)
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
color.tp53[color.tp53 == "FALSE"] <- "#75507BFF"
color.tp53[color.tp53 == "TRUE"] <- "#73D216FF"
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
color.prc <- c( "TRUE"="#73D216FF","FALSE" = "#75507BFF","grey93")
# annot.column$PRC2.inact
#Tissuelvs cell line
color.tissuevscl <- c("Tumor"="#D15A80", CL="#B9EDDC")

#Benign vs Malignant
color.bm <- c( Benign = "#4C9B6B", Malignant = "#8A2450")

#NF1 inact
color.nf1 <- as.character(annot.column$NF1.inact)
color.nf1[is.na(color.nf1)] <- "grey93"
color.nf1[color.nf1 == "FALSE"] <- "#75507BFF"
color.nf1[color.nf1 == "TRUE"] <- "#73D216FF"
color.nf1 <- unique(color.nf1)
if(length(color.nf1)==3){
  names(color.nf1) <- c("TRUE","FALSE","")
}else{
  names(color.nf1) <- c("TRUE","FALSE")
  
}


#CDKN2A Hotspot
color.cdkn2aH <- as.character(annot.column$CDKN2AH.inact)
color.cdkn2aH[is.na(color.cdkn2aH)] <- "grey93"
color.cdkn2aH[color.cdkn2aH == "FALSE"] <- "#75507BFF"
color.cdkn2aH[color.cdkn2aH == "TRUE"] <- "#315A09"
color.cdkn2aH <- unique(color.cdkn2aH)
# names(color.cdkn2aH) <- c("FALSE","","TRUE")
if(length(color.cdkn2aH)==3){
  names(color.cdkn2aH) <- c("FALSE","","TRUE")
}else{
  names(color.cdkn2aH) <- c("FALSE","TRUE")
  
}

color.cdkn2a <- as.character(annot.column$CDKN2A.inact)
color.cdkn2a[is.na(color.cdkn2a)] <- "grey93"
color.cdkn2a[color.cdkn2a == "FALSE"] <- "#75507BFF"
color.cdkn2a[color.cdkn2a == "TRUE"] <- "#73D216FF"
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
mycolors <- list(Sample.type = sample.name.col,
                 SpvsNF1 = color.nf1sp,
                 Seq.Platform = col.paltform,
                 TumorvsCL= color.tissuevscl,
                 Tumor.purity = color.purity,
                 BenignvsMalign = color.bm,
                 NF1 =color.nf1,
                 CDKN2A =color.cdkn2a,
                 CDKN2A.SV.hotspot =color.cdkn2aH,
                 PRC2 = color.prc,
                 TSG_inactivation = color.tsg,
                 TP53 = color.tp53,
                 Activating.mut = color.actmut,
                 CN.LOH.0.5 = color.cn.loh,
                 LOH.0.5 = color.wgs.loh,
                 Mut.freq = color.mut.freq,
                 Even.chr.50 = color.even.chr,
                 Methylome_ID = color.methyl,
                 Genomic.group.order = col.groups,
                 H3K27me3 = col.H3k27me3,
                 Clusters = col.clusters,
                 chr1.unbalance = color.chr1.unb
                 
)
pvalue <- 0.05

annot.column[sapply(annot.column, is.logical)] <- lapply(annot.column[sapply(annot.column,is.logical)],as.character)
annot.column[sapply(annot.column, is.character)][sapply(annot.column[sapply(annot.column, is.character)],is.na)] <-""

#Transcription factors genes
tf.genes <-read.table("TFsGenes_Lovering2020.txt")[,1]

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
rdata.file.rlog <- file.path("./DataFigures","dds.rlog_AllSamplesAllGenes_RealSet.RData")

if(file.exists(rdata.file.rlog)){
  load(file.path("./DataFigures","dds.rlog_AllSamplesAllGenes_RealSet.RData"))
}else{
  dds.rlog <- rlog(dds.deg)
  save(dds.rlog, file = file.path("./DataFigures","dds.rlog_AllSamplesAllGenes_RealSet.RData"))
  
}


#dds.rlog TF 
dds.rlogs <- dds.rlog[rownames(dds.rlog)%in%tf.genes,]


################ PCA representation ###################

pca.pseudocounts <- dds.rlogs
#Intgroups are the colnames present in colData(dds.rlogs)
intgroups <- c("Sample.type","SpvsNF1","TumorvsCL","BenignvsMalign","Tumor.purity","NF1","CDKN2A","TP53","PRC2","Clusters","Genomic.group.order")
i=11
for(i in seq_len(length(intgroups))){
  intgroup <- intgroups[i]
  if(intgroup%in%c("NF1","CDKN2A","TP53","PRC2","Clusters")){
    colData(pca.pseudocounts)[c(1:9,30:32,41:49),intgroup]<-NA
    
  }
  
  if(intgroup == "CDKN2A"){
    sv.hotspot <- "CDKN2A.SV.hotspot"
    sv.hotspot <- sample.data[,sv.hotspot]
    names(sv.hotspot) <- sample.data$Sample.name
    sv.hotspot <-sv.hotspot[names(sv.hotspot)%in% rownames(colData(pca.pseudocounts))]
    sv.hotspot <- data.frame(sv.hotspot)
    sv.hotspot$sv.hotspot[is.na(sv.hotspot$sv.hotspot)] <-FALSE
    # sv.hotspot[c(1:9,30:32,41:49),"sv.hotspot"]<-NA
    sv.hotspot.true <- rownames(sv.hotspot)[sv.hotspot$sv.hotspot==TRUE]
    colData(pca.pseudocounts)[sv.hotspot.true,intgroup] <- "SV"
    mycolors$CDKN2A <- c(mycolors$CDKN2A,"SV" ="#315A09")
  }
  pcaData <- plotPCA(pca.pseudocounts, intgroup=intgroup, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  color.sample <- mycolors[[intgroup]]
 
  rownames(pcaData)
  if(intgroup %in% c("Tumor.purity")){
    pal.col <- brewer.pal(n=9,"Blues")
    # brewer.blues
    # brewer.pal.info
    # png(filename=file.path(pca.dir, paste0("PCA_PC1-PC2_TF_AllSamples_GenesPurity.IdentitySamples",sample.group,"nonames_IntBiol_RealSet.png")), width = 2000, height = 1800)
    # svg(filename=file.path(pca.dir, paste0("PCA_PC1-PC2_TF_AllSamples_GenesPurity.IdentitySamples",sample.group,"nonames_IntBiol_RealSetValsetTesis_v3.svg")), width = 40, height = 28)
    
    pp <- ggplot(pcaData, aes(PC1, PC2, color = pcaData[,intgroup]))+
      geom_point(size = 45, stroke =3)+ 
      scale_color_gradient(name = "%Tumor Cells", low ="#DEEBF7", high =  "#08306B")+
      # ggrepel::geom_label_repel(aes(label = as.character(rownames(pcaData))),size =10, show.legend = FALSE)+
      theme(legend.title = element_text(size=60),
            legend.key.size =unit(3,"cm"),
            legend.position = "top",
            legend.text = element_text( size=40,), 
            axis.text=element_text(size=70),
            axis.title=element_text(size=70),
            plot.title = element_text(size=50),
            panel.background = element_rect(fill="white"), 
            panel.grid.major = element_line(colour = "grey50"), plot.margin =margin(30,100,30,30))+
      guides(shape = guide_legend( override.aes = list(size = 30, color ="#7C8e8F", fill = "#7C8e8F")))+
      xlab(paste0("PC1")) +
      ylab(paste0("PC2"))
    ggsave(plot = pp,filename =file.path(pca.dir, paste0("PCA_PC1-PC2_TF_AllSamples_Genes.Purity.IdentitySamples",sample.group,"nonames_IntBiol_RealSet_v4.svg")), width =40, height = 28)
    
  }else{
    pp <- ggplot(pcaData, aes(PC1, PC2, color = pcaData[,intgroup] ))+
      geom_point(size = 45, stroke =3,pch = 16)+ 
      scale_color_manual(name = intgroup, values =color.sample)+
      # ggrepel::geom_label_repel(aes(label = as.character(rownames(pcaData))), show.legend = FALSE, size = 10)+
      theme(legend.title = element_text(size=60), 
            legend.text = element_text( size=60), 
            legend.position = "top",
            axis.text=element_text(size=70),
            axis.title=element_text(size=70),
            plot.title = element_text(size=50),
            panel.background = element_rect(fill="white"), 
            panel.grid.major = element_line(colour = "grey50"), plot.margin =margin(30,100,30,30))+
      guides(colour = guide_legend(override.aes = list(shape = 16, size = 15)),
             shape = guide_legend( override.aes = list(size = 5, color ="#7C8e8F", fill = "#7C8e8F" )))+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance"))
    ggsave(plot = pp,filename =file.path(pca.dir, paste0("PCA_PC1-PC2_TFgenes_AllSamples.", intgroup,".", sample.group,".nonames_IntBiol_RealSet_v4.svg")), width =40, height = 28)
    # ggsave(plot = pp,filename =file.path(pca.dir, paste0("PCA_PC2-PC3_SelectedGenesItzi_OnlyTumorSamples.", intgroup,".", sample.group,".nonames_IntBiol_RealSet_v4.svg")), width =40, height = 28)
    
  }
}





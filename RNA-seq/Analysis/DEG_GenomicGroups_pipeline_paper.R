#####################################################
#    Integrative biology: DEG genomic groups        #
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
file.names <- sample.data$Paper.names
sample.names <- sample.data$Sample.name
sample.type <-unique(sample.data[,sample.group][sample.data[,sample.group]!=""])

#Salmon parameters
output.suffix <- "_quant"
output.quants <- "./results/Salmon"
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


#Transcription factors genes
tf.genes <-read.table("./DataFigures/TFsGenes_Lovering2020.txt")[,1]


# table(geness %in% genes)
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

#Saving txi.salmon$abundance (TPM) for OTTER classifier
tpm.df <-data.frame(txi.salmon$abundance)
tpm.df$gene_name <- rownames(tpm.df)
tpm.df <- tpm.df[,c(ncol(tpm.df),1:(ncol(tpm.df)-1))]
colnames(tpm.df)[2:ncol(tpm.df)] <- colnames(txi.salmon$abundance)

samples <- colnames(tpm.df)
samples <- samples[-1]
for(i in seq_len(length(samples))){
  sn <- samples[i]
  df_sn <- data.frame(gene_name = rownames(tpm.df),
                      TPM = tpm.df[,sn])
  write.table(df_sn, file = file.path("results/Revision", paste0(sn,"_TPM_table_OTTER.tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
}



write.table(tpm.df, file = file.path("results/Revision/TPM_table_OTTER.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

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
}else{
  dds.rlog <- rlog(dds.deg)
  save(dds.rlog, file = rdata.file.rlog)
  
}
colnames(dds.rlog)

#dds.rlog TF 
dds.rlogs <- dds.rlog[rownames(dds.rlog)%in%genes,]
dds.rlogs <- dds.rlogs[,colData(dds.rlogs)$Genomic.group.order %in% c(1,2,3,4)]
#################      Differential gene expression      ########################
#Sample.data 
sample.sel.data <- sample.data[sample.data[,sample.group]!="", ]
sample.sel.data <- sample.sel.data[grepl("Miller",sample.sel.data$Sample.name),]
sample.sel.data <- sample.sel.data[sample.sel.data$Genomic.group.order %in% c(1,2,3,4),]

sample.sel.data$Genomic.group.order[sample.sel.data$Genomic.group.order%in%c(1,2)] <- "G1G2"
sample.sel.data$Genomic.group.order[sample.sel.data$Genomic.group.order%in%c(3,4)] <- "G3G4"
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
sample.type <- unique(sample.sel.data[,sample.group])
# sample.type <- sample.type[-c(5,6)]

#Getting the Deseq results
msg(paste("Starting DESeq of", sample.group, "contrast"))

# deseq.results.file <- file.path(DEG.dir,"DESeqResults_Genomic.group.order_G1G2vsG3G4_Valset_AllGenes.RData")
deseq.results.file <- file.path(DEG.dir,"DESeqResults_Genomic.group.order_Valset_AllGenes.RData")

if(file.exists(deseq.results.file)){
  load(deseq.results.file)
  
}else{
  deseq.results <- list()
  for(i in seq_len(length(sample.type)-1)) {
    # for(i in seq_len(length(con1))) {
    ref <- as.character(sample.type)[i]
    # ref <- con1[i]
    
    deseq.results[[ref]] <- list()
    
    for(j in (i+1):length(sample.type)) {
      # for(j in 1:length(con2)) {
      other_group <- as.character(sample.type)[j]
      # other_group <- con2[j]
      
      # FIPS samples and time %in% c(ref, otro)
      t <- c(ref, other_group)
      deseq_samp <- sample.sel.data[sample.sel.data[,sample.group] %in% t, c("Sample.name", sample.group)]
      
      #runDESeq function
      
      message("RunDeSEQ ", other_group, " vs ", ref)
      
      deseq.results[[ref]][[other_group]] <-runDESeq( sample_names = deseq_samp$Sample.name, sname_variable = "Sample.name", ref_group = ref, samples_df = sample.sel.data, filtered_dds = filtered.dds, type_lfcShrink = "apeglm", test = "Wald")
      # deseq.results[[ref]][[other_group]] <-res
      msg("runDESeq done")
      
    }
    
  }
  
  save(deseq.results,file = deseq.results.file)
  
}

# Getting Differentially expressed genes

selected_genes <- list()
tosave <-F
for( i in 1:length(deseq.results)){
  lists <- deseq.results[[i]]
  class(lists)
  ref <- names(deseq.results[i])
  contrast <- data.frame()
  for (j in 1:length(lists)){
    contrast <- lists[[j]]
    other_group <- names(lists[j])
    
    # We select those genes which have a pvale <= 0.05 from all analised genes.
    selected_genes[[ref]][[other_group]]<- getTopGenesLFC(deseq_results = contrast, padj_value = 0.05)
    
    selected_genes[[ref]][[other_group]]$genes <- rownames(selected_genes[[ref]][[other_group]])
    # selected_genes[[ref]][[other_group]] <- selected_genes[[ref]][[other_group]][rownames(selected_genes[[ref]][[other_group]]) %in% tf.genes,]
    
    # We save the results obtained
    if(tosave){
      write.table(selected_genes[[ref]][[other_group]],
                  file.path(selected.genes.dir, paste0("ResultTable_", ref,"vs", other_group,"_",sample.group,"_ValSet.csv")),
                  sep = "\t", dec = ".", quote = FALSE, col.names = TRUE, row.names = FALSE)
      
      
      write.table(selected_genes[[ref]][[other_group]]$genes,
                  file.path(selected.genes.dir, paste0("DEGenes_", ref,"vs", other_group,".DEG.",sample.group,"_ValSet.csv")),
                  sep = "\t", col.names = FALSE, row.names = FALSE)
      
    }
   
  }
}

genes <- selected_genes$`1`$`3`$genes
################################# Enrichment DEG ##############################
genes <- selected_genes$G1G2$G3G4$genes[selected_genes$G1G2$G3G4$log2FoldChange <(-1)]
lfc.g <- selected_genes$G1G2$G3G4$log2FoldChange[selected_genes$G1G2$G3G4$log2FoldChange >(0)]
names(lfc.g) <- genes
lfc.genes <- sort(lfc.g,decreasing = TRUE)

genes <- selected_genes$`1`$`2`$genes[selected_genes$`1`$`2`$log2FoldChange>=(1)]
lfc.g <- selected_genes$`1`$`2`$log2FoldChange[selected_genes$`1`$`2`$log2FoldChange >=(1)]
names(lfc.g) <- genes
lfc.genes <- sort(lfc.g,decreasing = TRUE)

genes <- selected_genes$`1`$`3`$genes[selected_genes$`1`$`3`$log2FoldChange>(0)]
genes <- genes[genes%in% selected_genes$`2`$`3`$genes[selected_genes$`2`$`3`$log2FoldChange>(0)]]


enrich.terms <- getAllEnrichTermList(names(lfc.genes))
enrich.terms$KEGG
for(i in seq_len(length(goterms))){
  t <- goterms[i]
  msg(t)
  term.path <-file.path(gsea.dir, t)
  if(nrow(data.frame(enrich.terms[[t]]))<1)next
  write.table(data.frame(enrich.terms[[t]]), file.path(term.path, paste0("Enrichment_", t, "_G1G2vsG3G4_G3G4lfc0_",sample.group,".v1.csv")),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  #### Enrichment plots
  barplot(enrich.terms[[t]],showCategory = 30, font.size = 4)
  ggsave(filename = file.path(term.path,paste0("Barplot_Enrich_",t,"_G1G2vsG3G4_G3G4lfc0.png")), width = 1000, height = 1000,units =  "px")
  
  clusterProfiler::cnetplot(x = enrich.terms[[t]], showCategory = 30, node_label = "gene", foldChange = lfc.genes, color.params= list(foldChange = lfc.genes), cex.params = list(gene_label=1.5,category_label = 0.5))+
    # clusterProfiler::cnetplot(x = enrich.terms[[i]][[t]], showCategory = 20,cex.params = list(gene_label=0.5,category_label = 0.5))+
    theme(panel.background = element_rect(fill = "white", colour = "white"), plot.background = element_rect(fill = "white"))
  # ggsave(filename = file.path(term.path,t,paste0("Enrich_",t,"_cluster_",ii,"_DEG_3KO_NCvs2KO_NC_cluster_",k,"_cnetplot_node.gene_noWT.v6.png")), width = 3500, height = 2000,units =  "px")
  ggsave(filename = file.path(term.path,paste0("Enrich_",t,"_G1G2vsG3G4_G3G4lfc0_cnetplot_node.gene.svg")), width = 15, height = 10)
  
  clusterProfiler::cnetplot(x = enrich.terms[[t]], showCategory = 30, node_label = "all", foldChange = lfc.genes, color.params= list(foldChange = lfc.genes),cex.params = list(gene_label=0.5,category_label = 1))+
    # clusterProfiler::cnetplot(x = enrich.terms[[i]][[t]], showCategory = 20,cex.params = list(gene_label=0.5,category_label = 0.5))+
    theme(panel.background = element_rect(fill = "white", colour = "white"), plot.background = element_rect(fill = "white"))
  ggsave(filename = file.path(term.path, paste0("Enrich_",t,"_G1G2vsG3G4_G3G4lfc0_cnetplot_node.all.svg")), width = 15, height = 10)
  
  clusterProfiler::cnetplot(x = enrich.terms[[t]], showCategory = 30, node_label = "category", foldChange = lfc.genes, color.params= list(foldChange = lfc.genes),cex.params = list(gene_label=1,category_label = 1))+
    # clusterProfiler::cnetplot(x = enrich.terms[[i]][[t]], showCategory = 20,cex.params = list(gene_label=0.5,category_label = 0.5))+
    theme(panel.background = element_rect(fill = "white", colour = "white"), plot.background = element_rect(fill = "white"))
  ggsave(filename = file.path(term.path, paste0("Enrich_",t,"_G1G2vsG3G4_G3G4lfc0_cnetplot_node.term.svg")), width = 15, height = 10)
  
 
  
}

 


################################# Representation DEG ##########################

##### Heatmap representation ######

#Selection of genes to represent
data_subset_norm <- dds.rlog[,colData(dds.rlog)$Genomic.group.order %in%c(1:4)]
data_subset_norm <- dds.rlogsel[,colData(dds.rlogsel)$Genomic.group.order %in%c("G1G2","G3G4")]

data_subset_norm <- assay(data_subset_norm)[rownames(data_subset_norm) %in%genes,]

#Column order
anc <- annot.column[rownames(annot.column) %in% colnames(data_subset_norm),] 
anc$Genomic.group <- factor(anc$Genomic.group,levels =c(1,2,3,4))
anc <- anc[order(anc$Genomic.group),]
sn.order <- rownames(anc)
data_subset_norm <- data_subset_norm[,sn.order] 

# Column annotation 
colnames(annot.column)
anno.col <- annot.column[,c(1,6,8,10,11,12,20,2)] 

png(filename = file.path(heatmap.dir,paste0("DEG_G1G2vsG3G4_upG1G2",sample.group,"_ColumnUNSupervised_IntBiol_Valset.v2.png")), width = 1200, height = 1000)
# png(filename = file.path(heatmap.dir,paste0("30PC3genes_supervised.2_AllSamples.IntBiol.png")), width = 1500, height = 2500)
# png(filename = file.path(heatmap.dir,paste0("MesenchymalMarkers_Itzi_",sample.group,"_UnSupervised.AllSamples.Realset.png")), width = 1200, height = 1000)

pheatmap(data_subset_norm,
         # annotation_row = anno.row,
         annotation_col = anno.col,
         annotation_colors = mycolors,
         cluster_cols = T,
         cluster_rows = T,
         # breaks = myBreaks,
         fontsize = 10, fontsize_row = 8,
         color = bluered(80),
         show_rownames =F,
         annotation_legend = TRUE, 
         legend = TRUE,
         annotation_names_row = F,
         annotation_names_col = T,
         treeheight_row = 50, 
         treeheight_col = 50,
         scale ="row",
         lwd = 1)
dev.off() 

##### Boxplot and barplot representation ####

genomic.groups.data <-sample.data[sample.data$Genomic.group.order %in% c(1,2,3,4),]
#Only Miller samples
genomic.groups.data <- genomic.groups.data[grepl("Miller",genomic.groups.data$Sample.name),]
txi.sel <- selectDataFromTximport(tximport = txi.salmon,sample_names = genomic.groups.data$Sample.name )
filtered.dds.groups <- getFilteredDDS(tximport = txi.sel, samples_group = sample.group, samples_df = genomic.groups.data, filter_min_reads = filt.min.reads, filter_min_samples = filt.min.samples)

filt.counts <- counts(filtered.dds.groups, normalize = TRUE)

#Selecting genes
filt.counts <- filt.counts[rownames(filt.counts) %in% c(genes),]
fct <- t(filt.counts)
fct <- data.frame(fct)

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

#Boxplot representation
for(i in seq_len(ncol(fct)-1)){
  # png(filename = file.path("/imppc/labs/eslab/mmagallon/Projects/cNTAP/cNTAP_DEG/Results/LuminexGenes/GenesLuminexCounts",
  # paste0("Boxplot_", colnames(fct)[i], "_normalizedCounts_ALL_luminexGenes_GroupSamples.png")),width = 1200, height = 1000)
  # 
  m <-max(max(fct[,i]))
  
  if(m<=999){
    plus <- 50
    
  } else{
    plus <- 500
    
  }
  pos <- round_any(seq(0,m+plus,by = plus), plus, f = ceiling)
  pos.at  <- c(min(pos),mean(pos),max(pos))
  
  png(filename = file.path(markers.dir,
                           paste0("Boxplot_", colnames(fct)[i], "_NOnormalizedCounts_GenomicGroupSamples_G1G2vsG3g4_G1G2_lfcup1.png")),width = 1200, height = 1000)
  
  par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")
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
                xaxt="n",yaxt="n", ylim= c(0, max(m)+plus ))
  
  axis(1, at = c(0:15),labels = rep("",16),lwd.ticks = 0, lwd = 3)
  
  
  # if(m>1000) d <--3
  # if(m<1000) d <- -2
  
  
  # pos <- round(seq(0, m+max(fct[,i]) , by = round(m/3,d=0)),d =0)
  # pos <- round(c(min(pos), mean(pos), max(pos)))
  # pos <- round_any(c(0, m/2, m), 500, f = ceiling)
  
  
  axis(2, at = pos.at, las=2, lwd = 3, cex.axis=5, ylim =c(0,max(m)))
  # axis(2, at = seq(0,m+100,50000),las=2, lwd = 3, cex.axis=5)
  # axis(2, at = pos, las = 2, lwd = 3, cex.axis =5,ylim =c(0,max(pos)) )
  legend("topright", inset=c(-0.4,0),box.lwd = "", xpd = TRUE, fill = unique(color.data), legend = levels(fct$stage), horiz = F, cex = 3,col = 2)
  # legend(x = 8,y=1000,box.lwd = "", xpd = TRUE, fill = color.data[colnames(matrix.data)], legend = names(color.data[colnames(matrix.data)]), horiz = F, cex = 1,col = 2)
  
  dev.off()
}

#Barplot representation

filt.counts <- filt.counts[,rownames(fct)]
i=8
for( i in seq_len(nrow(filt.counts))){
  png(filename = file.path(file.path(markers.dir, paste0("Barplot_",rownames(filt.counts)[i],"_GenomicGroupSamples_G1G2vsG3g4_G1G2_lfcup1.png"))),
      width = 1500, height = 1000)
  barplot(filt.counts[i,],las=2,col = color.data, main= rownames(filt.counts)[i])
  dev.off()
}



df <- data.frame(fct)
colnames(df)
dft <- df[df$stages %in% c(2,1),]
genes.df <- colnames(df)
genes.df <- genes.df[-length(genes.df)]
g.pv <- c()
for(i in seq_len(length(genes.df))){
  g <- genes.df[i]
  tt <- t.test(dft[,g]~dft$stages)
  g.pv <- c(g.pv, tt$p.value)
  names(g.pv)[i] <-g
}
g.pv.02 <- data.frame("pval"=g.pv[g.pv<0.05])
enrichment <- getAllEnrichTermList(rownames(g.pv.02))
write.table(g.pv.02, file = file.path("./results/Genomic.group.order/markers/Mirar/C6", "genes_pval_0.2.csv"),sep = "\t",col.names = T, row.names = T)
lm1 <- lm(df$ERG8~stages,data = df)
anno <- aov(lm1)
TukeyHSD(anno)

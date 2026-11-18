#########################################
#     PCA projection Validation Set     #
#########################################
#Packages needed
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(gplots)
library(karyoploteR)

#Loading Data to plot
#### loading the sample data information
sample.data <- read.table(file = "./Sample.info.csv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

sample.data <- sample.data[sample.data$Real.ValSet %in% c("Real","ValSet"),]
sample.data$Genomic.group.order <- as.character(sample.data$Genomic.group.order)
sample.data$Genomic.group.order[is.na(sample.data$Genomic.group.order )] <- ""
sample.group <- "Clusters"
pca.dir <- file.path("results",sample.group,"PCAplot")
################## Color of samples ################
annot.column <- data.frame(sample.name = sample.data$Sample.type,
                           SpvsNF1 =sample.data$SpvsNF1,
                           TumorvsCL = sample.data$TumorvsCL,
                           BeningvsMalign = sample.data$BenignvsMalign,
                           # Purity = sample.data$Tumor.purity,
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
                           # Mut.freq = sample.data$Mut.freq,
                           Even.chr = sample.data$Even.chr.50,
                           # platform = sample.data$Seq.Platform,
                           Methylation.id =sample.data$Methylome_ID,
                           Genomic.group = sample.data$Genomic.group.order,
                           Cluster = sample.data$Clusters,
                           Sample.set = sample.data$Real.ValSet,
                           chr1.unbalance = sample.data$chr1.unbalance)

rownames(annot.column) <- sample.data$Sample.name
#Samplet type colors
sample.name.col <- unique(sample.data$Group.color)
names(sample.name.col) <- unique(sample.data$Sample.type)
sample.name.col <- unique(sample.name.col)
names(sample.name.col) <- unique(sample.data$Sample.type)

#Platform colors
col.paltform <- c(ILL = "yellow", BGI ="orange")

# #Purity color ramp
# purity <- annot.column$Purity
# pal <- colorRampPalette(c("darkblue","white"))
# color.purity <- map2color(x = purity,pal =pal(15) ,limits = c(0,1))
# names(color.purity) <-unique(sample.data$Sample.name)
# color.purity[is.na(color.purity)] <- "grey"
# 
# #Mut.freq
# mut.freq <- annot.column$Mut.freq
# pal <- colorRampPalette(c("white","orange"))
# color.mut.freq <- map2color(x = mut.freq,pal =pal(19) ,limits = c(0,3e5))
# names(color.mut.freq) <-unique(sample.data$Sample.name)
# color.mut.freq[is.na(color.mut.freq)] <- "grey"

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
color.tp53[color.tp53 == "FALSE"] <-  "#75507BFF"
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
  names(color.nf1) <- c("FALSE","TRUE")
  
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
color.chr1.unb <- c( "TRUE"="#47B4F7","FALSE" = "grey72", "grey93")
#87CEFA
#clusters
col.clusters <- c(C1= "#7BA7FFFF", C2 = "#FFEC6DFF","Grey"= "grey","grey92")
#col.set
col.set <- c(Real = "#0B007A",ValSet= "#87CEFA")

annot.column$LOH.prop[annot.column$LOH.prop == ""] <- "FALSE"
annot.column$Even.chr[annot.column$Even.chr == ""] <- "TRUE"
mycolors <- list(Sample.type = sample.name.col,
                 Subtype = color.nf1sp,
                 # platform = col.paltform,
                 TumorvsCL= color.tissuevscl,
                 # Purity = color.purity,
                 BeningvsMalign = color.bm,
                 NF1 =color.nf1,
                 CDKN2A =color.cdkn2a,
                 CDKN2A.SV.hotspot =color.cdkn2aH,
                 PRC2 = color.prc,
                 TSG_inactivation = color.tsg,
                 TP53 = color.tp53,
                 Activating.mut = color.actmut,
                 CN.LOH.0.5 = color.cn.loh,
                 LOH.0.5 = color.wgs.loh,
                 # Mut.freq = color.mut.freq,
                 Even.chr = color.even.chr,
                 Methylation.id = color.methyl,
                 Genomic.group.order = col.groups,
                 H3K27me3 = col.H3k27me3,
                 Clusters = col.clusters,
                 Real.ValSet = col.set,
                 chr1.unbalance = color.chr1.unb
                 
)
pvalue <- 0.05

annot.column[sapply(annot.column, is.logical)] <- lapply(annot.column[sapply(annot.column,is.logical)],as.character)
annot.column[sapply(annot.column, is.character)][sapply(annot.column[sapply(annot.column, is.character)],is.na)] <-""



# Data to plot in PCA projection
tf.genes <-read.table("TFsGenes_Lovering2020.txt")[,1]

load(file.path("./DataFigures","dds.rlog_AllSamplesAllGenes_RealSet.RData"))
dds.rlog.real <- dds.rlog[rownames(dds.rlog)%in%tf.genes,]
pcaData.real <- plotPCA(dds.rlog.real, intgroup=c("Sample.type"), returnData=TRUE)

load(file.path("./DataFigures","dds.rlog_AllSamplesAllGenes_RealSetValSet.RData"))
dds.rlog.val <- dds.rlog[rownames(dds.rlog)%in%tf.genes,]

dds.rlog.val <- dds.rlog.val[,50:99]
colnames((dds.rlog.val))
# assay(dds.rlog.val)[, colnames(assay(ddsgrepl(colData(dds.rlog.val),"Miller"))]


rv <- rowVars(assay(dds.rlog.real))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,
                                                   length(rv)))]


intgroups <- c("Sample.type","PRC2","NF1","TSG_inactivation","Subtype","Clusters",
               "Genomic.group.order","H3K27me3","chr1.unbalance","CDKN2A","CDKN2A.SV.hotspot","Activating.mut","CN.LOH.0.5", "Real.ValSet")
# 
# 
# intgroup=c("PRC2")
# # intgroup=c("Methylome_ID")
# intgroup=c("Sample.type")
# # intgroup = "Main.Dx.MPNST"
# # intgroup = "Histological.features"
# intgroup = "TSG_inactivation"
# intgroup = "Subtype"
# intgroup = "Scoring"
# intgroup = "Clusters"
# intgroup = "Genomic.group.order"
# intgroup = "H3K27me3"
# intgroup = "chr1.unbalance"
# intgroup = "CDKN2A.SV.hotspot"
# intgroup = "Activating.mut"
# intgroup ="Mut.freq"
# intgroup = "CN.LOH.0.5"
# sample.data$
i=6
for(i in seq_len(length(intgroups))){
  intgroup <- intgroups[i]
  pca = prcomp(t(assay(dds.rlog.real)[select, ]))
  project.b = predict(pca,t(assay(dds.rlog.val)))
  pcax <- pca$x
  pcaData <- rbind(pcax,project.b)
  
  intgroup.df <- as.data.frame(colData(dds.rlog.real)[, intgroup,
                                                      drop = FALSE])
  intgroup.df.2 <- as.data.frame(colData(dds.rlog.val)[, intgroup,
                                                       drop = FALSE])
  intgroup.df <- rbind(intgroup.df,intgroup.df.2)
  intgroup.df[c(1:9,30:32,41:49),intgroup]<-NA
  if(intgroup == "CDKN2A"){
    sv.hotspot <- "CDKN2A.SV.hotspot"
    sv.hotspot <- sample.data[,sv.hotspot]
    names(sv.hotspot) <- sample.data$Sample.name
    sv.hotspot <-sv.hotspot[names(sv.hotspot)%in% rownames(intgroup.df)]
    sv.hotspot <- data.frame(sv.hotspot)
    sv.hotspot$sv.hotspot[is.na(sv.hotspot$sv.hotspot)] <-FALSE
    # sv.hotspot[c(1:9,30:32,41:49),"sv.hotspot"]<-NA
    sv.hotspot.true <- rownames(sv.hotspot)[sv.hotspot$sv.hotspot==TRUE]
    intgroup.df[sv.hotspot.true,intgroup] <- "SV"
    mycolors$CDKN2A <- c(mycolors$CDKN2A,"SV" ="#315A09")
  }
  pcaData <- cbind(pcaData, 
                   intgroup.df, name = rownames(intgroup.df))
  percentVar <- round(100*(pca$sdev^2/sum(pca$sdev^2)))
  color.sample <- mycolors[[intgroup]]
  class(color.sample)
  
  
  # svg(filename=file.path(pca.dir, paste0("PCA_PC1-PC2_TFgenes_AllSamples.", intgroup,".", sample.group,".nonames_IntBiol_RealSetValSet_projection_v4.svg")), width = 40, height = 28)
  pp <- ggplot(pcaData, aes(PC1, PC2, color = pcaData[,intgroup] ))+
    geom_point(size = 45, stroke =3,pch = 16)+ 
    scale_color_manual(name = intgroup, values =color.sample)+
    # ggrepel::geom_label_repel(aes(label = as.character(rownames(pcaData))), show.legend = FALSE)+
    theme(legend.title = element_text(size=60), 
          legend.text = element_text( size=60), 
          legend.position = "top",
          axis.text=element_text(size=70),
          axis.title=element_text(size=70),
          plot.title = element_text(size=50),
          panel.background = element_rect(fill="white"), 
          panel.grid.major = element_line(colour = "grey50"), plot.margin =ggplot2::margin(30,100,30,30))+
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 15)),
           shape = guide_legend( override.aes = list(size = 5, color ="#7C8e8F", fill = "#7C8e8F" )))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  ggsave(plot = pp,filename =file.path(pca.dir, paste0("PCA_PC1-PC2_TFgenes_AllSamples.", intgroup,".", sample.group,".nonames_IntBiol_RealSetValSet_projection_v5.svg")), width =40, height = 28  )
  # dev.off()
  j=1
  if(intgroup == "Genomic.group.order"){
    for(j in seq_len(length(names(col.groups)))){
      gr <- names(col.groups)[j]
      mycolors$Genomic.group <- col.groups
      mycolors$Genomic.group[length( mycolors$Genomic.group)] <- "grey"
      # mycolors$Genomic.group <- mycolors$Genomic.group[-(length(mycolors$Genomic.group)-1)]
      mycolors$Genomic.group[!names(mycolors$Genomic.group)%in%c(gr,"")] <- transparent("grey72")
      pcaData$Genomic.group.order <- as.character(pcaData$Genomic.group.order)
      pcaData$Genomic.group.order[is.na(pcaData$Genomic.group.order)]<-""
      
      # svg(filename=file.path(pca.dir, paste0("PCA_PC1-PC2_TFgenes_AllSamples.GenomicsGroup_",gr,"_",sample.group,"_nonames_IntBiol_RealSetValSet_projection_v3.svg")), width = 40, height = 28)
      
      pp <- ggplot(pcaData, aes(PC1, PC2, color = Genomic.group.order ))+
        geom_point(size = 45, stroke =3,pch = 16)+ 
        scale_color_manual(name = "Genomic group", values =mycolors$Genomic.group)+
        # ggrepel::geom_label_repel(aes(label = as.character(rownames(pcaData))), show.legend = FALSE)+
        theme(legend.title = element_text(size=60),
              legend.text = element_text( size=60),
              legend.position = "top",
              axis.text=element_text(size=70),
              axis.title=element_text(size=70),
              plot.title = element_text(size=50),
              panel.background = element_rect(fill="white"),
              panel.grid.major = element_line(colour = "grey50"), plot.margin =ggplot2::margin(30,100,30,30))+
        guides(colour = guide_legend(override.aes = list(shape = 16, size = 15)),
               shape = guide_legend( override.aes = list(size = 5, color ="#7C8e8F", fill = "#7C8e8F" )))+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      ggsave(plot = pp, filename =file.path(pca.dir, paste0("PCA_PC1-PC2_TFgenes_AllSamples.",intgroup,"_",gr,"_",sample.group,"_nonames_IntBiol_RealSetValSet_projection_v5.svg")), width =40, height = 28  )
      
    }
  }
}




library(karyoploteR)
library(CopyNumberPlots)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library("BSgenome.Hsapiens.UCSC.hg38.masked")
# library(reticulate)

#Function
plotCopyNumberCalls <- function (karyoplot, cn.calls, cn.values = NULL, cn.column = "cn", 
                                 cn.colors = NULL, loh.values = NULL, loh.column = "loh", 
                                 loh.color = "#1E90FF", loh.height = 0.3, labels = NULL, label.cex = 1, 
                                 label2.cex = NULL, track.margin = 0.01, r0 = 0, r1 = 1, ...) 
{
  if (methods::is(cn.calls, "list") || methods::is(cn.calls, 
                                                   "CompressedGRangesList")) {
    labels <- prepareLabels(labels = labels, x = cn.calls)
    for (i in seq_len(length(cn.calls))) {
      rr <- karyoploteR::autotrack(current.track = i, total.tracks = length(cn.calls), 
                                   margin = track.margin, r0 = r0, r1 = r1)
      plotCopyNumberCalls(karyoplot, cn.calls[[i]], r0 = rr$r0, 
                          r1 = rr$r1, cn.values = cn.values, loh.values = loh.values, 
                          cn.column = cn.column, loh.column = loh.column, 
                          cn.colors = cn.colors, loh.color = loh.color, 
                          labels = labels[i], label.cex = label.cex, label2.cex = label2.cex, 
                          loh.height = loh.height, ...)
    }
    return(invisible(karyoplot))
  }
  if (is.null(label2.cex)) 
    label2.cex <- label.cex
  cn.calls <- regioneR::toGRanges(cn.calls)
  if (is.null(labels)) 
    labels <- "CN"
  if (!is.null(labels) && !is.na(labels) && all(is.character(labels) && 
                                                length(labels) > 0)) {
    if (length(labels) == 1) {
      karyoploteR::kpAddLabels(karyoplot, labels = labels[1], 
                               r0 = r0, r1 = r1, cex = label.cex, ...)
    }
    else {
      karyoploteR::kpAddLabels(karyoplot, labels = labels[1], 
                               r0 = r0 + (r1 - r0) * loh.height, r1 = r1, cex = label.cex, 
                               ...)
      karyoploteR::kpAddLabels(karyoplot, labels = labels[2], 
                               r0 = r0, r1 = r0 + (r1 - r0) * loh.height, cex = label2.cex, 
                               ...)
    }
  }
  if (is.null(cn.values)) {
    if (length(GenomicRanges::mcols(cn.calls)) == 0) 
      stop("No cn.values given and cn.calls has no associated copy number data")
    if (is.null(cn.column)) {
      cn.values <- names(GenomicRanges::mcols(cn.calls))[1]
    }
    else {
      if (!(cn.column %in% names(GenomicRanges::mcols(cn.calls)))) 
        stop("The cn.calls object does not have a column ", 
             cn.column, ". No copy number data is available")
      cn.values <- GenomicRanges::mcols(cn.calls)[, cn.column]
    }
  }
  if (is.null(loh.values)) {
    if (loh.column %in% names(GenomicRanges::mcols(cn.calls))) {
      loh.values <- GenomicRanges::mcols(cn.calls)[, loh.column]
    }
  }
  if (!is.null(loh.values)) {
    if (!is.logical(loh.values)) 
      loh.values <- tryCatch(as.logical(loh.values))
    loh.values[is.na(loh.values)] <- FALSE
  }
  segment.colors <- getCopyNumberColors(cn.colors)
  plot.region <- karyoplot$plot.region
  segments <- IRanges::subsetByOverlaps(cn.calls, plot.region)
  seg.cols <- segment.colors[as.character(cn.values)]
  karyoploteR::kpRect(karyoplot, data = cn.calls, y0 = loh.height, 
                      y1 = 1, col = seg.cols, r0 = r0, r1 = r1, border = NA, 
                      ...)
  if (!is.null(loh.values)) {
    karyoploteR::kpRect(karyoplot, data = cn.calls[loh.values], 
                        y0 = 0, y1 = loh.height, r0 = r0, r1 = r1, col = loh.color, 
                        border = NA, ...)
  }
  invisible(karyoplot)
}

########## Execution #####


execution.dir <- "./execution.Dir"

sample.data <- read.table(file = file.path(execution.dir, "Sample.info.txt"), sep ="\t", header = T, comment.char = "")
sample.names <- sample.data$Sample.name
# sample.names <- sample.names[c(1,2,3,6,7,14:36)]
# sample.names <- sample.names[grepl("BLOOD", sample.names)]
# sample.data$Purity.CNVkit[33]<-0.6
# sample.data$Purity.CNVkit[32]<- 0.85
# sample.data$Purity.CNVkit[36]<- 0.4
# sample.data$Ploidy
genome <- "hg38"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
orgdb <- org.Hs.eg.db


load("./CN_AllMPNSTcellLinespaper_LOHcorrected_annotated.v5.RData")


sample.cn <- sample.cn[2:4]
names(sample.cn)
# genes <- c("CDKN2C","CKS1B","FGFR3","CHIC2","PDGFRA","NUP98","KMT2A","ATM","TP53")
genes <- c("CDKN2C","CKS1B","NUP98","KMT2A")

gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

tx <- GenomicFeatures::genes(x= TxDb.Hsapiens.UCSC.hg38.knownGene)
entrezid2symbol <- biomaRt::select(org.Hs.eg.db,mcols(tx)$gene_id, c("SYMBOL"))
tx$SYMBOL <- entrezid2symbol$SYMBOL
gene.markers.gf <- tx[tx$SYMBOL%in%genes]
gene.markers.gf <- sort(gene.markers.gf,ignore.strand = TRUE)

############## Figures 
pp <- getDefaultPlotParams(plot.type = 4)
pp$topmargin <- 10
pp$bottommargin <- 20
pp$data1inmargin <- 10
pp$ideogramheight <- 3
pp$leftmargin <- 0.20
pp$rightmargin <- 0.05
pp$data1outmargin <- 30
pp$ideogramlateralmargin <- 0.001

plotDefaultPlotParams(plot.type = 3,pp)


cn.colors <- getCopyNumberColors()
cn.colors[3] <- "#F0F0F0"
# cn.colors <- transparent(cn.colors)

loh.calls <- unlist(as(sample.cn,"GRangesList"))
# loh.calls <- unlist(as(loh.calls,"GRangesList"))

lapply(sample.cn,is.null)
names(sample.cn)
which(is.null(sample.cn))
lapply(sample.cn, length)
loh.calls$loh <- as.logical(loh.calls$loh)
loh.calls <-loh.calls[loh.calls$loh ==TRUE] 
loh.calls <- loh.calls[,c("loh","cn","Sample.name")] 
chr <- "canonical"
label.axis <- c("0","50%","100%")
ymax = 1
cex.axis = 1.8
cex.titel.axis=2.5

#GENES to PLOT
i=1
# chrs <- c(seqlevels(sample.cn$`Tumor-18`))
# chrs <- unique(seqnames(essential.genes.reg))
chrs <- "canonical"
chr <- "canonical"
reg.def <- toGRanges(c("chr1:792501-121605200","chr4:69201-190048600","chr9:670035-39500291","chr10:45701-132557735","chr11:194501-135075600",
                       "chr14:80144882-106883718","chr16:46423001-90100500","chr17:404091-35766416"))

svg(filename=file.path(execution.dir,"results","images", paste0(chr,"_CNVprofile_MPNSTCellLines.svg")), width = 35, height = 15)
kp <- karyoploteR::plotKaryotype(genome = "hg38", plot.type = 4, chromosomes = chr, plot.params = pp, ideogram.plotter = NULL, labels.plotter = NULL)
kpPlotRegions(kp,data = reg.def,r0=0,r1=0.8, col = "#CCCDD4",border = NA)
karyoploteR::kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, cex = 2.5, srt = 45 )
kpPlotMarkers(kp, data = gene.markers.gf, labels = gene.markers.gf$SYMBOL, cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 0.8,y = 1, marker.parts = c(0.98, 0.05, 0.002),ignore.chromosome.ends = TRUE,clipping = FALSE)

ats <- autotrack(current.track = 1, total.tracks = 1, margin = 0.25)
plotCopyNumberCalls(kp, cn.calls = sample.cn , loh.height = 0.1, loh.column = "loh",  label.cex = 2.5 ,r0 =ats$r0, r1=ats$r1, cn.colors = cn.colors,  label.margin = 0.05)

# ats <- autotrack(current.track = 2, total.tracks = 2, margin = 0.25)
# 
# at <- autotrack(current.track = 1,total.tracks = 3, margin = 0.15,r0 = ats$r0,r1=ats$r1)
# plotCopyNumberSummary(kp, cn.calls = sample.cn , label.cex = 2.5,labels = "",r0 =at$r1, r1=at$r0,
#                       normal.color = "white", gain.color = NA, loss.color = "darkgreen", data.panel = 1)
# karyoploteR::kpAddLabels(karyoplot = kp,labels = "Losses", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
# 
# kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 =at$r1, r1=at$r0)
# kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
# 
# at <- autotrack(current.track = 2,total.tracks = 3, margin = 0.15,r0 = ats$r0,r1=ats$r1)
# plotCopyNumberSummary(kp, cn.calls = sample.cn , label.cex = 2.5,labels = "",r0 =at$r1, r1=at$r0,
#                       normal.color = "white",loss.color = NA, gain.color ="#FFD700",data.panel = 1,)
# karyoploteR::kpAddLabels(karyoplot = kp,labels = "Gains", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
# 
# kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = c("100%", "50%","0"), cex = cex.axis,r0 =at$r1, r1=at$r0)
# kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
# 
# at <- autotrack(current.track = 3,total.tracks = 3, margin = 0.15,r0 = ats$r0,r1=ats$r1)
# kpPlotCoverage(kp, data = loh.calls,
#                col = "#1C8FFFFF", data.panel = 1,
#                r0 =at$r0, r1=at$r1, ymax = length(sample.cn))
# 
# kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 = at$r0, r1=at$r1)
# kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 = at$r0, r1=at$r1)
# karyoploteR::kpAddLabels(karyoplot = kp,labels = "LOH", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
dev.off()

#FISH probes

# fish.probes <- list("chr1"= c("p32.3", "q21"),"chr4"=c("p16.3","q12"),
#                  "chr11"=c("p15.4","q23.3","q22.3"),"chr17"=c("p13","p11.1"))

fish.probes <- list("chr1"= c("p32.3", "q21"),
                    "chr11"=c("p15.4","q23.3"))


# chr.band <- lapply(sample.cn, function(x){
#   chr.band <- lapply(strsplit(",",x=x$chr.band),function(y) unique(y))
#   x$chr.band <- unlist(lapply(chr.band,function(y)paste0(y,collapse = ",")))
#            }
#   )
# for(i in seq_len(length(sample.cn))){
#   sn <- names(sample.cn)[i]
#   sample.cn[[sn]]$chr.band <- chr.band[[sn]]
# }



pp <- getDefaultPlotParams(plot.type = 2)
pp$topmargin <- 10
pp$bottommargin <- 20
pp$data1inmargin <- 10
pp$ideogramheight <- 30
pp$leftmargin <- 0.20
pp$rightmargin <- 0.05
pp$data1outmargin <- 30
pp$data1height <- 5
pp$ideogramlateralmargin <- 0.008

plotDefaultPlotParams(plot.type = 2,pp)


cn.colors <- getCopyNumberColors()
cn.colors[3] <- "#F0F0F0"

cytobands <- karyoploteR::getCytobands(genome = "hg38")
cytobands.gr <- c()
i=1
for(i in seq_len(length(fish.probes))){
  chr <- names(fish.probes)[i]
  bands <- fish.probes[[chr]]
  cytobands.chr <- cytobands[seqnames(cytobands)==chr]
  
  for(b in seq_len(length(bands))){
    bnd <- bands[b]
    
    cytobands.gr <- c(cytobands.gr,cytobands.chr[grepl(bnd,cytobands.chr$name)])
    
  }
}
cytobands.gr <- sort(unlist(as(cytobands.gr,"GRangesList")))
cytobands.gr <- sort(unlist(as(cytobands.gr,"GRangesList")))
cytobands.gr <- reduce(cytobands.gr)
cytobands.gr <- sort(toGRanges(cytobands.gr,genome = "hg38"))
cytobands.gr$cytobands <- unlist(fish.probes)
kp$plot.region
svg(filename=file.path(execution.dir,"results","images", paste0("Cytobands_MPNSTCellLines_chr11.svg")), width = 35, height = 15)
kp <- karyoploteR::plotKaryotype(genome = "hg38",plot.type = 2,chromosomes = "chr11",plot.params = pp,cex =6)
kpPoints(kp,data = c(gene.markers.gf),y = -8,cex = 8,col = transparent(c("green","green")))
kpPoints(kp,chr = "chr11",x = end(gene.markers.gf)+1000000,y = -8,cex = 8,col = transparent(c("red","red")))

# kpText(kp,data=cytobands.gr,y = -3,labels = cytobands.gr$cytobands, cex =4)
kpPlotMarkers(kp, data = gene.markers.gf, labels = gene.markers.gf$SYMBOL, cex =4, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.05, 0.002),ignore.chromosome.ends = TRUE,clipping = FALSE)
plotCopyNumberCalls(kp, cn.calls = sample.cn,  cn.colors = cn.colors,label.cex = 4,data.pane =2)
# kpPoints()
dev.off()

s <-1
for(s in seq_len(length(cytobands.gr))){
  band <- cytobands.gr[s]
  # sn <- names(sample.cn)[s]
  # sn.cn <- sample.cn[[sn]]
  
  svg(filename=file.path(execution.dir,"results","images", paste0(band,"_Cytobands_MPNSTCellLines.svg")), width = 15, height = 12)
  kp <- karyoploteR::plotKaryotype(genome = "hg38",plot.type = 4,zoom = band,plot.params = pp, main = paste0(seqnames(band),band$cytobands),cex = 2)
  kpPlotMarkers(kp, data = gene.markers.gf, labels = gene.markers.gf$SYMBOL, cex =1.2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.05, 0.002),ignore.chromosome.ends = TRUE,clipping = FALSE)
  
  plotCopyNumberCalls(kp, cn.calls = sample.cn, cn.colors = cn.colors,label.cex = 1.5)
  dev.off()
  
}




    ploidy.file <- file.path(execution.dir, "results", sample, "CNVkit/ExcludedRegions2_ploidy/cnvkit_test_ploidy", paste0(sample,"_call_ploidy_",pd,"_purity_", pur,".cns"))


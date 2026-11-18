######################
# Figure CDKN2A ANFs #
######################
#Pakages needed
library(karyoploteR)
library(CopyNumberPlots)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library("BSgenome.Hsapiens.UCSC.hg38.masked")
# library(reticulate)
library(biomaRt)
biomart.dataset <- useEnsembl(biomart= "genes", dataset="hsapiens_gene_ensembl",version = 108)

##### Functions
#

# plotCopyNumberCalls <- function (karyoplot, cn.calls, cn.values = NULL, cn.column = "cn", 
#                                  cn.colors = NULL, loh.values = NULL, loh.column = "loh", 
#                                  loh.color = "#1E90FF", loh.height = 0.3, labels = NULL, label.cex = 1, 
#                                  label2.cex = NULL, track.margin = 0.01, r0 = 0, r1 = 1, ...) 
# {
#   if (methods::is(cn.calls, "list") || methods::is(cn.calls, 
#                                                    "CompressedGRangesList")) {
#     labels <- prepareLabels(labels = labels, x = cn.calls)
#     for (i in seq_len(length(cn.calls))) {
#       rr <- karyoploteR::autotrack(current.track = i, total.tracks = length(cn.calls), 
#                                    margin = track.margin, r0 = r0, r1 = r1)
#       plotCopyNumberCalls(karyoplot, cn.calls[[i]], r0 = rr$r0, 
#                           r1 = rr$r1, cn.values = cn.values, loh.values = loh.values, 
#                           cn.column = cn.column, loh.column = loh.column, 
#                           cn.colors = cn.colors, loh.color = loh.color, 
#                           labels = labels[i], label.cex = label.cex, label2.cex = label2.cex, 
#                           loh.height = loh.height, ...)
#     }
#     return(invisible(karyoplot))
#   }
#   if (is.null(label2.cex)) 
#     label2.cex <- label.cex
#   cn.calls <- regioneR::toGRanges(cn.calls)
#   if (is.null(labels)) 
#     labels <- "CN"
#   if (!is.null(labels) && !is.na(labels) && all(is.character(labels) && 
#                                                 length(labels) > 0)) {
#     if (length(labels) == 1) {
#       karyoploteR::kpAddLabels(karyoplot, labels = labels[1], 
#                                r0 = r0, r1 = r1, cex = label.cex, ...)
#     }
#     else {
#       karyoploteR::kpAddLabels(karyoplot, labels = labels[1], 
#                                r0 = r0 + (r1 - r0) * loh.height, r1 = r1, cex = label.cex, 
#                                ...)
#       karyoploteR::kpAddLabels(karyoplot, labels = labels[2], 
#                                r0 = r0, r1 = r0 + (r1 - r0) * loh.height, cex = label2.cex, 
#                                ...)
#     }
#   }
#   if (is.null(cn.values)) {
#     if (length(GenomicRanges::mcols(cn.calls)) == 0) 
#       stop("No cn.values given and cn.calls has no associated copy number data")
#     if (is.null(cn.column)) {
#       cn.values <- names(GenomicRanges::mcols(cn.calls))[1]
#     }
#     else {
#       if (!(cn.column %in% names(GenomicRanges::mcols(cn.calls)))) 
#         stop("The cn.calls object does not have a column ", 
#              cn.column, ". No copy number data is available")
#       cn.values <- GenomicRanges::mcols(cn.calls)[, cn.column]
#     }
#   }
#   if (is.null(loh.values)) {
#     if (loh.column %in% names(GenomicRanges::mcols(cn.calls))) {
#       loh.values <- GenomicRanges::mcols(cn.calls)[, loh.column]
#     }
#   }
#   if (!is.null(loh.values)) {
#     if (!is.logical(loh.values)) 
#       loh.values <- tryCatch(as.logical(loh.values))
#     loh.values[is.na(loh.values)] <- FALSE
#   }
#   segment.colors <- getCopyNumberColors(cn.colors)
#   plot.region <- karyoplot$plot.region
#   segments <- IRanges::subsetByOverlaps(cn.calls, plot.region)
#   seg.cols <- segment.colors[as.character(cn.values)]
#   karyoploteR::kpRect(karyoplot, data = cn.calls, y0 = loh.height, 
#                       y1 = 1, col = seg.cols, r0 = r0, r1 = r1, border = NA, 
#                       ...)
#   if (!is.null(loh.values)) {
#     karyoploteR::kpRect(karyoplot, data = cn.calls[loh.values], 
#                         y0 = 0, y1 = loh.height, r0 = r0, r1 = r1, col = loh.color, 
#                         border = NA, ...)
#   }
#   invisible(karyoplot)
# }

########## Parameters ##########
execution.dir <- "executionDir"

sample.data <- read.table(file = file.path(execution.dir, "Sample.info.csv"), sep ="\t", header = T,comment.char = "")
sample.names <- sample.data$Sample.name[-c(1,30:32)]
if(is.character(sample.data$Purity.CNVkit)){
  sample.data$Purity.CNVkit <- gsub(",",".",sample.data$Purity.CNVkit)
  sample.data$Ploidy <- gsub(",",".",sample.data$Ploidy)
}
sample.names <- sample.names[sample.names %in% sample.data$Sample.name[sample.data$Cluster == "aNF"]]
names(sample.names) <- c("ANNUBP-1","","ANNUBP-4","ANNUBP-5","ANNUBP-2","ANNUBP-6","ANNUBP-7","ANNUBP-3","","","")
sample.names <- sample.names[names(sample.names) != ""]
genome <- "hg38"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
orgdb <- org.Hs.eg.db

ref.genome.fq <- "./GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref.genome <- "hg38"
ploidy <- 2

#Loading CN variants
sample.log <- list()
sample.cn <- list()
i=2
for(i in seq_len(length(sample.names))){
  sample <- sample.names[i]
  sn <- sample.names[i]
  #obtainin the purity of de sample previously calculated
  pur <- sample.data$Tumor.purity.est[sample.data$Sample.name == sample]
  # pd <- sample.data$Ploidy[sample.data$Sample.name == sample]
  # pur <- purity
  if(!sn %in% names(sample.log)){ 
    # Data to plot
    log.file <- file.path(execution.dir,"results", sample,"CNVkit/ExcludedRegions/TumorvsPooledNormals", paste0(sample,".cnr"))
    logs <-loadCopyNumberCallsCNVkit(log.file)
    logs <- filterChromosomes(logs)
    names(mcols(logs))[2] <- "lrr"
    sample.log[[sn]] <-logs
  }
  
  
  for(k in seq_len(length(ploidy))){
    pd <- ploidy[k]
    # sample.cn[[sn]][[paste0("ploidy_",pd)]] <- loadCopyNumberCallsCNVkit(file.path(execution.dir, "results", sample, "CNVkit/ExcludedRegions/cnvkit_test_ploidy", paste0(sample,"_call_ploidy_",pd,".cns")))
    ploidy.file <- file.path(execution.dir, "results", sample, "CNVkit/ExcludedRegions/TumorvsPooledNormals/cnvkit_test_ploidy", paste0(sample,"_call_ploidy_",pd,"_purity_", pur,".cns"))
    # ploidy.file <- file.path(execution.dir, "results", sample, "CNVkit/ExcludedRegions/TumorvsPooledNormals/cnvkit_test_ploidy", paste0(sample,"_call_ploidy","_purity_", pd,".cns"))
    if(!file.exists(ploidy.file))next
    
    sample.cn[[sn]][[paste0("ploidy_",pd)]] <- loadCopyNumberCallsCNVkit(ploidy.file)
    sample.cn[[sn]][[paste0("ploidy_",pd)]] <- filterChromosomes(sample.cn[[sn]][[paste0("ploidy_",pd)]])
    
    # loh
    diff <- which(sample.cn[[sn]][[paste0("ploidy_",pd)]]$cn1 != sample.cn[[sn]][[paste0("ploidy_",pd)]]$cn2)
    cn1 <- sample.cn[[sn]][[paste0("ploidy_",pd)]]$cn1[diff]
    cn2 <- sample.cn[[sn]][[paste0("ploidy_",pd)]]$cn2[diff]
    sel <- cn1 !=0 & cn2 == 0
    sel.2 <- cn1 ==0 & cn2!=0
    diff <- c(diff[sel],diff[sel.2])
    
    sample.cn[[sn]][[paste0("ploidy_",pd)]]$loh <-FALSE
    sample.cn[[sn]][[paste0("ploidy_",pd)]]$loh[diff] <- TRUE
    
  }
  
  
}

lumpy.file <- "All_Lumpy_SV_Regions_ANFs_FilteredCommonSVinSamples_mw6.v7.RData"
load(file.path(execution.dir,"results/", lumpy.file))
names(filt.list)



#########################################################
#             Chromosomal representation                #
#########################################################
# Plot parameters
cn.colors <- getCopyNumberColors()
cn.colors[3] <- "white"

pp <- getDefaultPlotParams(plot.type = 1)
pp$topmargin <- 30
pp$bottommargin <- 20
pp$data1inmargin <- 5
pp$ideogramheight <- 2
pp$ideogramlateralmargin <- 0.003
pp$leftmargin <- 0.1
plotDefaultPlotParams(plot.type =1, plot.params = pp)

gene.markers <- read.table(file.path("./","special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
gene.markers.gf <- toGRanges(gene.markers$V2, genome ="hg38")
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6

load(file.path(execution.dir,"results","ANFs_BAF_profile.RData"))
#### CDKN2A regions ####

cdkn2aARF.gene <- toGRanges("chr9:21967733-21999999")
cdkn2a.gene <- toGRanges("chr9:21967733-21975870")
locus.cdkn2a <- toGRanges("chr9:21779106-22202743")
cdkn2a.3mb <- toGRanges("chr9:18967753-24975098")
cdkn2a.2mb <- toGRanges("chr9:19967753-23975098")
cdkn2a.1mb <- toGRanges("chr9:21550000-22437372")
regions <- c(cdkn2a.gene=cdkn2a.gene, cdkn2aARF.gene=cdkn2aARF.gene, locus.cdkn2a=locus.cdkn2a, cdkn2a.3mb=cdkn2a.3mb, cdkn2a.2mb=cdkn2a.2mb, cdkn2a.1mb=cdkn2a.1mb)

# Getting cdkn2a region
kp <- plotKaryotype(zoom = toGRanges("chr9:21.98e6-21.99e6"), genome="hg38")

gg <- karyoploteR::makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, karyoplot = kp, plot.transcripts.structure = TRUE)
gg <- addGeneNames(gg)
gg$genes <- gg$genes[gg$genes$name=="CDKN2A"] #Remove MTAP, CDKN2B-AS, CDKN2A-DT
cdkn2a.tr <- gg$transcripts$`1029`
cdkn2a.tr$tx_id
data.frame(cdkn2a.tr)
# cdkn2a.tr <- cdkn2a.tr[cdkn2a.tr$tx_id %in% c(105800, 105801)] #Keep only the two main transcripts
# cdkn2a.tr <- cdkn2a.tr[cdkn2a.tr$tx_id %in% c(114866, 114867)] #Keep only the two main transcripts
cdkn2a.tr <- cdkn2a.tr[cdkn2a.tr$tx_name %in% c("ENST00000304494.10",	"ENST00000579755.2")] #Keep only the two main transcripts

#cdkn2a.tr <- cdkn2a.tr[!(cdkn2a.tr$tx_id %in% c("105812", "105811", "105809", "105803", "105813"))]
gg$transcripts$`1029` <- cdkn2a.tr

tt <- list("1029"=gg$transcripts$`1029`)
gg$transcripts <- tt


trans  <- cdkn2a.tr$tx_id
ce <-list()
nce <- list()
i=1
for(i in seq_len(length(trans))){
  tr <- as.character(trans[i])
  ce[[tr]] <- gg$coding.exons[[tr]]
  nce[[tr]] <- gg$non.coding.exons[[tr]]
  
}
gg$coding.exons <- ce
gg$non.coding.exons <- nce

CDKN2A.reg <- gg$genes[gg$genes$name=="CDKN2A"]
plot.zoom <- CDKN2A.reg+1e4



kp <- plotKaryotype(genome="hg38", zoom=plot.zoom, ideogram.plotter = NULL, plot.params = pp, cex=2, labels.plotter = NULL)
kpAddChromosomeNames(kp, cex=5,xoffset = -0.05)
#kpAddCytobandsAsLine(kp, lwd = 4)
kpAbline(kp, h=0, lwd=4, data.panel="ideogram", col="#444444")
kpAddBaseNumbers(kp, tick.dist = 15000, minor.tick.dist = 2000, add.units = TRUE, cex =8, tick.len = 6, minor.tick.len = 3, lwd=2)



#gg.merged
kpPlotGenes(kp, data = gg , r1=0.5, add.gene.names = FALSE, add.transcript.names = TRUE, transcript.names = c("105800"=expression(p16^INK4a), "105801"=expression(p14^ARF)), transcript.name.position = "left", transcript.name.cex = 4, 
            introns.col = "#555555", mark.distance = 6, mark.height = 0.2, lwd=5)
kpText(kp, chr="chr9", x=21.9625e6, y=0.2, labels = "CDKN2A", cex=8, clipping=FALSE)

#Add the exon names
kpText(kp, gg$coding.exons$`119157`[1], labels = "E3", y=-0.025, pos=1, clipping=FALSE, cex=8, offset=1)
kpText(kp, gg$coding.exons$`119157`[2], labels = "E2", y=-0.025, pos=1, clipping=FALSE, cex=8, offset=1)
kpText(kp, gg$coding.exons$`119157`[3], labels = "E1\u03B1", y=-0.025, pos=1, clipping=FALSE, cex=8, offset=1)
kpText(kp, gg$coding.exons$`119158`[2], labels = "E1\u03B2", y=-0.025, pos=1, clipping=FALSE, cex=8, offset=1)

filt.list$PNF1_22$SelectedRegionsFiltered
filt.list$`66ANF`$SelectedRegionsFiltered <- c(filt.list$`66ANF`$SelectedRegionsFiltered,toGRanges(data.frame(chr = "chr9", start=21971515, end =21971515, REF = "N", ALT = "[chr9:21971515[N", 
                                                                                                              IMPRECISE = "TRUE", SVTYPE = "ICR", SVLEN = NA, END = "chr20:30500964", MATEID = "<NA>", STRANDS = "", Sample= "66ANF", Pos.var= "chr9:21971515-21971515",
                                                                                                              repearMasker = "",SYMBOL = "CDKN2A", DEL.ext.genes = "",DUP.et.genes = "", Mate.Genes = "",X.hg38.= "hg38")))
                                               
filt.list$`66ANF`$SelectedRegionsFiltered <- c(filt.list$`66ANF`$SelectedRegionsFiltered,toGRanges(data.frame(chr = "chr9", start=21971676, end =21971676, REF = "N", ALT = "N]chr9:21971676]",
                                                                                                              IMPRECISE = "TRUE", SVTYPE = "ICR", SVLEN = NA, END = "chr20:33262003", MATEID = "<NA>", STRANDS = "", Sample= "66ANF", Pos.var= "chr9:21971676-21971676",
                                                                                                              repearMasker = "",SYMBOL = "CDKN2A", DEL.ext.genes = "",DUP.et.genes = "", Mate.Genes = "",X.hg38.= "hg38")))

                                                                                                              
filt.list$PNF4II_2$SelectedRegionsFiltered
#########
#CDKN2A_ARF
kp <- plotKaryotype(plot.type = 1,  plot.params = pp, zoom = regions[[2]],
                      cex= 3, genome = "hg38",ideogram.plotter = NULL, labels.plotter = NULL)
tick.dist <- 8e3

# #1e6
# kp <- plotKaryotype(plot.type = 1,  plot.params = pp, zoom = regions[[6]],
#                     cex= 3, genome = "hg38",ideogram.plotter = NULL, labels.plotter = NULL)
# tick.dist <- 2e5
# 
# kp <- plotKaryotype(plot.type = 1,  plot.params = pp, zoom = regions[[1]],
#                     cex= 3, genome = "hg38",ideogram.plotter = NULL, labels.plotter = NULL)
# tick.dist <- 1e3

clipping <- TRUE
add.units=FALSE
digits=2

toLabel <- function(n, add.units, digits) {
  if(add.units==TRUE) {
    units <- c("b", "Kb", "Mb")
  } else {
    units <- c("", "", "")
  }
  if(abs(n) < 1000) return(paste0(as.character(n), units[1]))
  if(abs(n) < 1000000) return(paste0(as.character(round(n/1000, digits=digits)), units[2])) #Kb
  return(paste0(as.character(round(n/1000000, digits=digits)), units[3])) #Mb
}

old.scipen <- options("scipen")
options(scipen=999)
on.exit(options(scipen=old.scipen), add=TRUE)


#For every chromsome

karyoplot <- kp
for(chr.name in karyoplot$chromosomes) {
  
  chr <- karyoplot$genome[chr.name]
  #Major ticks
  num.ticks <- width(chr)/tick.dist + 1
  
  tick.pos <- start(chr) + (tick.dist*(0:(num.ticks-1))) - 1 
  tick.pos[1] <- start(chr)
  
  #if zoomed in, keep only the ticks in the plot region
  if(karyoplot$zoom==TRUE) {
    if(clipping==TRUE) {
      tick.pos <- tick.pos[tick.pos >= start(karyoplot$plot.region) & tick.pos<= end(karyoplot$plot.region)]
    }
  }
}


#1mbCDKN2A
# base.numbers <- toGRanges(data.frame(chr= "chr9", start=tick.pos, end = tick.pos, label = c(21.6,21.8,22,22.2,22.4)))

#CDKN2AARF
base.numbers <- toGRanges(data.frame(chr= "chr9", start=tick.pos, end =tick.pos, label = c(21.96,21.97,21.98,21.99)))


sample.names

windows <- 5e3
# max.dens <- 5447
# max.dens <- 1981
# max.dens <- 1.09

# max.dens <- 1981
#when we look at 1mb around CDKN2A we will select max dens of 500
max.dens <- 500
max.dens <- 80
sample.names <- sample.names[c(1,4,7,2,3,5,6)]
sample.names <- sample.names[length(sample.names):1]
r <- transparent("darkorchid")
i=2
# for(i in seq_len(length(regions))){
  region  <- regions[[i]]
  
  # png(filename = file.path("./MPNST_tumors/WGS/results/images/", paste0("CDKN2A_ANFs_SV_BAMlarge_tumorsGRAND_NAMES_PAPER.", names(regions[i]),"_max.dens_",max.dens,".png")), width = 2500, height =1200)
  # png(filename = file.path(images.dir.SV, paste0("CDKN2A_MPNST_SV_tumorsGRAND_NAMES_CHR9_max.dens_",max.dens,".png")), width = 3500, height = 2000)
  svg(filename = file.path("./MPNST_tumors/WGS/results/images/", paste0("CDKN2A_ANFs_SV_BAMlarge_tumorsGRAND_NAMES_PAPER.", names(regions[i]),"_max.dens_",max.dens,".svg")), width = 40, height =18)
  
  # region = "chr9"
  # Plotting CDKN2A locus
  kp <- plotKaryotype(plot.type = 1,  plot.params = pp, zoom = region, cex= 3, genome = "hg38",ideogram.plotter = NULL, labels.plotter = NULL)
  # s <-kpAddBaseNumbers(kp, 2e5, cex = 2, clipping = TRUE)
  # kpAddBaseNumbers(kp, 8e3, cex = 2, clipping = TRUE)
  if(names(regions)[i]== "cdkn2a.1mb"){
    kpPlotRegions(kp, data = toGRanges(c("chr9:21967733-21995870")), col = "#9932CC80", r0=0, r1 = 1.4)
    
  }
  
  # kp <- plotKaryotype(chromosomes = region, plot.type = 1,  plot.params = pp, cex= 4, genome = "hg38")
  t.tracks <- 11
  
  #Plot Genes
  at <-autotrack(current.track = 1:2, total.tracks = t.tracks)
  kpPlotGenes(kp, gg,  data.panel = 1, add.gene.names = TRUE,
              add.transcript.names = FALSE, r0=at$r0, r1=at$r1, 
              gene.name.cex =3, gene.name.position = "left", clipping = TRUE, mark.height = 0.5 )
  
  
  
  at<- autotrack(current.track = 3, total.tracks = t.tracks)
  kpSegments(kp,data = base.numbers, y0 =0.7, y1 =1, r0=at$r0, r1=at$r1,col = "black", lwd=2)
  kpText(kp, data = base.numbers, labels = base.numbers$label, pos =1, y = 0.7, r0=at$r0, r1=at$r1, clipping = FALSE,  cex =2.5)
  
  
  # Coverage 
  # max.density <- c()
  # windows <- 5e3
  windows <- 1e3
  # at <- autotrack(current.track = 2:t.tracks, total.tracks = t.tracks, margin = 0.1)
  j=1
  svs <- c()
  for(j in seq_len(length(sample.names))){
    sn <- sample.names[j]
    if(!sn%in% names(filt.list))next
    sn.sv <- filt.list[[sn]][["SelectedRegionsFiltered"]]
    del <- sn.sv[sn.sv$SVTYPE == "DEL"]
    transloc <-  sn.sv[sn.sv$SVTYPE %in% c("ICR")]
    svs <- c(svs,transloc)
    
    inv <-  sn.sv[sn.sv$SVTYPE %in% c("INV")]
    svs <- c(svs,inv)
    
    ats <-autotrack(current.track = t.tracks-j+1,total.tracks = t.tracks, margin = 0.25)
    plotCopyNumberCalls(kp, cn.calls = sample.cn[[sn]], r0=ats$r0, r1=ats$r1, loh.values = FALSE, loh.height = 0, labels ="", cn.colors = cn.colors)
    # kpPlotBAMDensity(kp, data = file.path("./MPNST_tumors/WGS/results",sn,"/BAM", paste0(sn,".bam")), r0=ats$r0, r1=ats$r1, window.size = windows, col = "grey", ymax = max.dens)
    kpPlotBAMCoverage(kp, data = file.path("./MPNST_tumors/WGS/results",sn,"/BAM", paste0(sn,".bam")), r0=ats$r0, r1=ats$r1, window.size = windows, col = "grey", ymax = max.dens)
    
    # if(sn=="38ANF1_old38PNF1")sn <- "38ANF1"
    # if(sn=="44ANFD_old44MPNSTD")sn <- "44ANFD"
    kpAddLabels(kp, labels = names(sn),r0 = ats$r0, r1 = ats$r1,cex = 3, label.margin = 0.03)
  
  }
  
  
  svs <- unlist(as(svs, "GRangesList"))
  kpPlotMarkers(kp, data = svs, labels = svs$END,
                line.color = "black" , cex =2.5, lty=2,
                lwd = 4, text.orientation = "vertical", label.margin = 2,
                srt = 45, label.dist = 0.03, max.iter = 1000, label.color = "black",#label.color = transloc$chr.col,
                r1 = 1.4, marker.parts = c(0.98,0.015,0.005), clipping = FALSE)
  
  for(j in seq_len(length(sample.names))){
    sn <- sample.names[j]
    if(!sn%in% names(filt.list))next
    sn.sv <- filt.list[[sn]][["SelectedRegionsFiltered"]]
    del <- sn.sv[sn.sv$SVTYPE == "DEL"]
    transloc <-  sn.sv[sn.sv$SVTYPE %in% c("ICR")]
    inv <-  sn.sv[sn.sv$SVTYPE %in% c("INV")]
    ats <-autotrack(current.track = t.tracks-j+1,total.tracks = t.tracks, margin = 0.25)
    
    if(length(transloc)>0){
      # kpText(kp, data= transloc, labels=transloc$END, y = 1.05,clipping = FALSE, pos =3, cex =2.5)
      # kpSegments(kp, data = transloc,  lty=2, y0 = ats$r0, y1 = 1.05, lwd = 4, clipping = FALSE)
      kpPoints(kp, data = transloc[transloc$Sample == sn], y = 0.1, pch = 17, r0=ats$r0, r1=ats$r1, col = "red" , cex = 3, lwd = 5)
      
    }
    if(length(inv)>0){
      # kpText(kp, data= transloc, labels=transloc$END, y = 1.05,clipping = FALSE, pos =3, cex =2.5)
      # kpSegments(kp, data = transloc,  lty=2, y0 = ats$r0, y1 = 1.05, lwd = 4, clipping = FALSE)
      kpPoints(kp, data = inv[inv$Sample == sn], y = 0.1, pch = 17, r0=ats$r0, r1=ats$r1, col = "blue" , cex = 3, lwd = 5)
      
    }
  }
  dev.off()
 
  # dev.off()
  
# }





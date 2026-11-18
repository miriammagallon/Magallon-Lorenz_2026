########################################################################## 
######## Integrative Biology of MPNST Tumors: CNV images  ################ 
########################################################################## 

##################### Packages needed ######################## 
library(CopyNumberPlots)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(biomaRt)
library(RColorBrewer)

######################## Functions ########################  
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


#Function kpPlotMeanCoveragebySample.
#This function was created for representing the mean of the number of alleles gained or lost
#We only add a parameter to kpPlotCoverage function from karyoploteR to easier represent the plot we wanted for
# cell lines characterization paper.

kpPlotMeanCoveragebySample <- function(karyoplot, data,
                                       num.samples, 
                                       show.0.cov=TRUE, 
                                       data.panel=1,
                                       r0=NULL, 
                                       r1=NULL, 
                                       col="#0e87eb", 
                                       border=NULL, 
                                       ymax=NULL, 
                                       clipping=TRUE, ...) {
  #Check parameters
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  #data
  if(missing(data)) stop("The parameter 'data' is required")
  #TODO: If data is not a SimpleRleList, try to convert it to a GRanges before testing with is so other Region Set formats can be used
  
  if(!methods::is(data, "GRanges") && !methods::is(data, "SimpleRleList")) {
    data <- tryCatch(toGRanges(data), error=function(e) {stop("'data' must be a GRanges object or a SimpleRleList")})
  }  
  
  #Compute (if needed) the coverage
  #If its not a coverage object,  it's a GRanges. Compute the coverage
  if(!methods::is(data, "SimpleRleList")) { 
    #remove any region not in the currently used genome
    data <- data[seqnames(data) %in% karyoplot$chromosomes,]
    #Old version, problems when data had no seqinfo - data <- GenomeInfoDb::keepSeqlevels(data, karyoplot$chromosomes, pruning.mode="coarse")
    #Remove any unused seq level from the GRanges to fix problems with coverage and witdh
    seqlevels(data) <- karyoplot$chromosomes
    #the width parameter is needed so the coverage extends to the end of the chromosomes
    data <- GenomicRanges::coverage(data, width=karyoplot$chromosome.lengths[seqlevels(data)]) 
  }
  
  coverage.gr <- toGRanges(data)
  
  if(show.0.cov==FALSE) {
    coverage.gr <- coverage.gr[coverage.gr$coverage!=0]
  }
  
  if(is.null(ymax)) ymax <- max(max(coverage.gr$coverage))
  
  if(is.null(border)) border <- col
  
  # kpBars(karyoplot=karyoplot, data=coverage.gr,
  #        y0=0, y1=coverage.gr$coverage, ymin=0, ymax=ymax,
  #        r0=r0, r1=r1, data.panel=data.panel,
  #        col=col, border=border, clipping=clipping, ...)
  
  
  #To get kpArea to plot the real coverage (flat tops), we need to build a
  #GRanges with two elements per range, one at the start and one at the end
  #NOTE: this breaks the show.cov.0=FALSE. 
  #TODO: Make kpArea and kpLines deal with NAs in the data and set coverage=0 
  #to NA here
  cov.start <- coverage.gr
  end(cov.start) <- start(cov.start)
  cov.end <- coverage.gr
  start(cov.end) <- end(cov.end)
  cov.to.plot <- sort(c(cov.start, cov.end))
  
  #mean of coverage gain by sample
  cov.to.plot$coverage <- cov.to.plot$coverage/num.samples
  max(cov.to.plot$coverage)
  kpArea(karyoplot=karyoplot, data=cov.to.plot, y=cov.to.plot$coverage, 
         base.y = 0, ymin=0, ymax=ymax,
         r0=r0, r1=r1, data.panel=data.panel,
         col=col, border=border, clipping=clipping, ...)
  
  karyoplot$latest.plot <- list(funct="kpPlotCoverage", computed.values=list(max.coverage=max(max(coverage.gr$coverage)),
                                                                             coverage=coverage.gr,
                                                                             ymax=ymax))
  
  invisible(karyoplot)
}

plotCopyNumberSummaryMeanCoverage <- function (karyoplot, cn.calls, direction = "in", gain.color = NULL, 
                                               normal.color = NULL, loss.color = NULL, add.grid = FALSE, 
                                               grid.color = "white", labels = NULL, label.cex = 1, label.srt = 0, 
                                               pos = 2, r0 = 0, r1 = 1, ...) {
  if (is.null(labels)) 
    labels <- "CN"
  if (!is.na(labels) && nchar(labels) > 0) {
    karyoploteR::kpAddLabels(karyoplot, labels = labels, 
                             r0 = r0, r1 = r1, cex = label.cex, srt = label.srt, 
                             pos = pos, ...)
  }
  direction <- match.arg(direction, c("in", "out"))
  if (any(is.null(gain.color), is.null(loss.color), is.null(normal.color))) {
    segment.colors <- getCopyNumberColors()
    if (is.null(gain.color)) 
      gain.color <- segment.colors["3"]
    if (is.null(normal.color)) 
      normal.color <- segment.colors["2"]
    if (is.null(loss.color)) 
      loss.color <- segment.colors["1"]
  }
  all.cn <- unlist(GenomicRanges::GRangesList(cn.calls))
  if (direction == "in") {
    karyoploteR::kpDataBackground(karyoplot, r0 = r0, r1 = r1, 
                                  col = normal.color)
    kpPlotMeanCoveragebySample(karyoplot, all.cn[all.cn$cn > 
                                                   2], num.samples = length(cn.calls), col = gain.color, r0 = r1, r1 = r0,  
                               ...)
    kpPlotMeanCoveragebySample(karyoplot, all.cn[all.cn$cn < 
                                                   2], num.samples = length(cn.calls), col = loss.color, r0 = r0, r1 = r1, 
                               ...)
    if (add.grid == TRUE) {
      karyoploteR::kpAbline(karyoplot, h = 0:length(cn.calls), 
                            col = grid.color, r0 = r0, r1 = r1, ymax = length(cn.calls), 
                            ...)
    }
  }
  else {
    mid <- r0 + (r1 - r0)/2
    kpPlotMeanCoveragebySample(karyoplot, all.cn[all.cn$cn > 
                                                   2], num.samples = length(cn.calls),col = gain.color, r0 = mid, r1 = r1, 
                               ...)
    kpPlotMeanCoveragebySample(karyoplot, all.cn[all.cn$cn < 
                                                   2], num.samples = length(cn.calls), col = loss.color, r0 = mid, r1 = r0, 
                               ...)
    if (add.grid == TRUE) {
      karyoploteR::kpAbline(karyoplot, h = 0:length(cn.calls), 
                            col = grid.color, r0 = mid, r1 = r1, ymax = length(cn.calls), 
                            ...)
      karyoploteR::kpAbline(karyoplot, h = 0:length(cn.calls), 
                            col = grid.color, r0 = mid, r1 = r0, ymax = length(cn.calls), 
                            ...)
    }
  }
  invisible(karyoplot)
}
setwd("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/")
source(file = "./Integrative_Biology_Functions.R")
source(file = "./CliffHunter.functions.R")
source(file ="/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/utils.R")
source(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-seq_MesenchimalDiff/RNA-seq_MesennchimalDiff/gsea_functions.R")

############### Parameters ########################
real.val.set <- "WGS"
# real.val.set <- "WGS_Validation"

execution.dir <- file.path("where/to/exevcute")
# setwd(execution.dir)

#ADD PARAMETER OF STRELKA DIRECTORY 
# strelka.dir <- input$strelka.dir
comm.res.dir <- "." 
if(real.val.set=="WGS"){
  # pseudobaf.rdata.file <- "AllIntegrativeTumors_WGS_BAF_Filtered.RData"
  pseudobaf.rdata.file <- "BAF_AllMPNSTandCellLinespaper_Filtered.RData"
  
  cnvkit.dir.raw <- "ExcludedRegions/TumorvsPooledNormals"
  sep = "\t"
  
  
}else if(real.val.set=="WGS_Validation"){
  #pseudobaf.rdata.file <- "TumorsMiler_BAF.RData"
  pseudobaf.rdata.file <- "BAF_AllMPNSTandCellLinespaper_Filtered.RData"
  sep = ","
  
}

sample.data <- read.table(file.path(execution.dir,"Sample.info.csv"), sep = sep, header = TRUE,comment.char = "")
if(real.val.set=="WGS"){
  sample.data <- sample.data[-c(1,30:32),]
  sample.names <- sample.data$Bioinfo.name
  names(sample.names) <- sample.data$Paper.Names
  sample.names <- sample.names[names(sample.names)!=""]
  sample.data <- sample.data[sample.data$Bioinfo.name %in% sample.names,]
}
if(real.val.set=="WGS_Validation"){
  sample.data <- sample.data[sample.data$BAM_Tumor ==TRUE,]
  sample.names <- sample.data$Sample.name
  names(sample.names) <- sample.data$Sample.name
  names(sample.names) <- gsub("Miller","ValSet",names(sample.names))
}


genome <- "hg38"

gene.markers <- toGRanges(data.frame(chr= c("chr9", "chr17", "chr17", "chr17", "chr11"), start = c(21967753,31094927, 31937007,7668230, 86245050), end = c(21975098, 31377677, 32001038, 7687366,86278810), symbol = c("CDKN2A", "NF1", "SUZ12", "TP53", "EED")),genome = genome)
info.param <- "avsnp150"
#Selection of interesting genes
gene.markers <- read.table(file.path(comm.res.dir,"special.genes.txt" ), 
                           header = FALSE, sep = " ", stringsAsFactors = FALSE) 
gene.markers.gf <- toGRanges(gene.markers$V2,genome = genome)
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6
regions.df <- toDataframe(regions.selected)
gene.markers.gf <- gene.markers.gf[gene.markers.gf$Genes !="BCR"]

load(file.path(execution.dir,"results", pseudobaf.rdata.file))
sns <- sample.data$Sample.name
names(sns)<- sample.data$Bioinfo.name
for(i in seq_len(length(sns))){
  sn <- sns[i]
  names(tumors.snp.lists)[names(tumors.snp.lists)%in%sn] <- names(sn)
}

# sample.names <- c(sample.names,"MPNST_SP_06" = "BN_MPNST")
if(is.character(sample.data$Tumor.purity.est)){
  sample.data$Tumor.purity.est <- gsub(",",".",sample.data$Tumor.purity.est)
  sample.data$Tumor.ploidy.est <- gsub(",",".",sample.data$Tumor.ploidy.est)
}
purity <- as.numeric(sample.data$Tumor.purity.est)
names(purity) <- names(sample.names)
class(purity)
purity <- sort(purity,FALSE)
col.purity <- colorRampPalette(brewer.pal(length(purity),"Blues"))(length(purity))
names(col.purity) <- names(purity)
sample.names <- sample.names[names(purity)]
col.purity <- col.purity[names(sample.names)]
sample.data$Sample.name 

#cluster to analize
c1 <- sample.data$Bioinfo.name[sample.data$Cluster %in% "C1"] 
c1 <- sample.names[sample.names %in% c1]
# c1 <- sample.names[names(sample.names)%in%c1]
# c1 <- c1[c1 != "MPNST-39_Miller"]
purity.c1 <- col.purity[names(c1)]
# c1[names(sample.names)%in%c1]
# c1 <- c1[length(c1):1] 
c2 <- sample.data$Bioinfo.name[sample.data$Cluster == "C2"] 
c2 <- sample.names[sample.names%in%c2]
purity.c2 <- col.purity[names(c2)]
sample.names <- c(c1,c2)
# sample.names <- sample.names[sample.names != "MPNST-39_Miller"]
# c2 <- c2[length(c2):1] 
# c2 <- c(c2,"MPNST_SP_06" = "BN_MPNST")

# CNVkit representation paper 

load(file.path("./CNV_data.RData"))
sample.cn <- sample.cn[sample.names] 
names(sample.cn) <- names(sample.names)
# names(sample.cn)==sample.names
i=1
for(i in  seq_len(length(sample.cn))){
  sn <- names(sample.cn)[i]
  sample.cn[[sn]]$Sample.name <-sn
}
# genes to plot 
#CNV UCSC track
all.genes.cosmic <- read.table(file.path("./Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]

tsg <- sort(all.genes.cosmic.gr[grepl("TSG", all.genes.cosmic$Role.in.Cancer)])
# seqlevels(tsg) <- seqlevels(sample.cn$`Tumor-18`)
tsg <- sort(tsg)

oncog <- sort(all.genes.cosmic.gr[grepl("oncogene", all.genes.cosmic$Role.in.Cancer)])
# seqlevels(oncog) <- seqlevels(sample.cn$`Tumor-18`)
oncog <- sort(oncog)

fusgenes <- sort(all.genes.cosmic.gr[grepl("fusion", all.genes.cosmic$Role.in.Cancer)])
# seqlevels(fusgenes) <- seqlevels(sample.cn$`Tumor-18`)
fusgenes <- sort(fusgenes)



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

########### Selection of SV regions >40% of samples in C1
loh.c1 <-loh.calls[loh.calls$Sample.name%in%c1]
loh.c1 <- loh.c1[loh.c1$cn == 2]
loh.c1 <- toGRanges(coverage(loh.c1))
max(loh.c1$coverage)
regions.genes <- loh.c1[loh.c1$coverage >5]  
reg.genes <- regions.genes[!seqnames(regions.genes)%in% c("chrX","chrY")]
reg.genes <- intersect(reg.genes,reg.genes)
reg.geness <- joinRegions(reg.genes,min.dist = 1e6)
reg.geness<- reg.geness[width(reg.geness)>1e6]
toDataframe(reg.genes)
reg.def <- toGRanges(c("chr1:792501-121605200","chr4:69201-190048600","chr9:670035-39500291","chr10:45701-132557735","chr11:194501-135075600",
                       "chr14:80144882-106883718","chr16:46423001-90100500","chr17:404091-35766416"))



loh.c2 <-loh.calls[loh.calls$Sample.name%in%c2]
loh.c2 <-subsetByOverlaps(loh.c2, reg.genes)


#annotation
biomart.dataset <- useEnsembl(biomart= "genes", dataset="hsapiens_gene_ensembl",version = 108)
reg.genes <- SVsBiomartGeneAnnotation(sv.data = reg.def, biomart.dataset)
reg.genes.df <- toDataframe(reg.genes)
reg.genes.df <- reg.genes.df[,-ncol(reg.genes.df)]

#Plot parameters

chr <- "chr9"
label.axis <- c("0","50%","100%")
ymax = 1
cex.axis = 1.8
cex.titel.axis=2.5

#GENES to PLOT
i=1
# chrs <- c(seqlevels(sample.cn$`Tumor-18`))
# chrs <- unique(seqnames(essential.genes.reg))
chrs <- "canonical"
i=1
for(i in seq_len(length(chrs))){
  chr <- chrs[i]
  res.dir <- file.path(execution.dir,"results","images")
  if(!file.exists(res.dir)) dir.create(res.dir)
  svg(filename=file.path(execution.dir,"results","images", paste0(chr,"_CNVprofile_cluster1.2_CN-LOHRegions_0.4_v6_Valset.svg")), width = 30, height = 20)
    total.tracks <- length(c1)+length(c2)+12
  kp <- karyoploteR::plotKaryotype(genome = "hg38", plot.type = 4, chromosomes = chr, plot.params = pp, ideogram.plotter = NULL, labels.plotter = NULL)
  
  # kpPlotMarkers(kp, data = oncog.loh, labels = oncog.loh$Gene.Symbol, cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.05, 0.002),ignore.chromosome.ends = TRUE,clipping = FALSE)
  # kpRect(kp,chr = seqnames(reg.genes), x0 = start(reg.genes),x1 = end(reg.genes),y0 = 0,y1=1,r0 = 0,r1=1, col = "#CCCDD4", border = NA,) 
  kpPlotRegions(kp,data = reg.def,r0=0,r1=1, col = "#CCCDD4",border = NA)
  karyoploteR::kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, cex = 2.5, srt = 45 )
  # at <- autotrack(current.track = 1:length(c2), total.tracks = total.tracks, margin = 0.20)
  # kpPlotMarkers(kp, data = tx.c3, labels = tx.c3$SYMBOL, cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.05, 0.002),ignore.chromosome.ends = TRUE,clipping = FALSE)
  
  # plotCopyNumberCalls(kp, cn.calls = sample.cn[names(c2)] , loh.height = 0.1, loh.column = "loh",  label.cex = 2.5 ,r0 =at$r0, r1=at$r1, cn.colors = cn.colors, track.margin = 0.25, label.margin = 0.05)
  for(r in seq_len(length(c2))){
    sn <- names(c2)[r]
    ats <- autotrack(current.track = r, total.tracks = total.tracks, margin = 0.25)
    plotCopyNumberCalls(kp, cn.calls = sample.cn[sn] , loh.height = 0.1, loh.column = "loh",  label.cex = 2.5 ,r0 =ats$r0, r1=ats$r1, cn.colors = cn.colors,  label.margin = 0.05)
    
    kpRect(kp,chr="chr1",x0=-150000000, x1 = -100000000,r0 = ats$r0, r1 = ats$r1, y0 =0, y1=1, clipping = TRUE, col =col.purity[names(c2)][r], track.margin = 0.25)
    # kpSegments(kp, chr = "chr1", x0=-150000000, x1 = 100000000,r0 =ats$r0, r1=ats$r1, y0 =0, y1=0)
    # kpSegments(kp, chr = "chr1", x0=-150000000, x1 = 100000000,r0 =ats$r0, r1=ats$r1, y0 =1, y1=1)
    
  }

  
  # rr <- karyoploteR::autotrack(current.track = 1:4, total.tracks = 6, 
  #                              r0 = at$r0, r1 = at$r1,margin = 0.25)
  # points(x = rep(0.15,length(names(c2))), y = seq(rr$r0+0.0,rr$r1+0.03,length.out=length(names(c2))), pch = 15,col = col.purity[names(c2)],cex =4)
  # 
  # abline(h =rr$r0+0.002)
  at <- autotrack(current.track = c(length(c2)+1,length(c2)+2),total.tracks = total.tracks, margin = 0.05)
  plotCopyNumberSummary(kp, cn.calls = sample.cn[names(c2)] , label.cex = 2.5,labels = "",r0 =at$r1, r1=at$r0,
                        normal.color = "white", gain.color = NA, loss.color = "darkgreen", data.panel = 1)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Losses", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 =at$r1, r1=at$r0)
  kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
  
  at <- autotrack(current.track = c(length(c2)+3,length(c2)+4),total.tracks = total.tracks, margin = 0.65)
  plotCopyNumberSummary(kp, cn.calls = sample.cn[names(c2)] , label.cex = 2.5,labels = "",r0 =at$r1, r1=at$r0,
                        normal.color = "white",loss.color = NA, gain.color ="#FFD700",data.panel = 1)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Gains", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = c("100%", "50%","0"), cex = cex.axis,r0 =at$r1, r1=at$r0)
  kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
  
  at <- autotrack(current.track = c(length(c2)+5,length(c2)+6),total.tracks = total.tracks, margin = 0.50)
  kpPlotCoverage(kp, data = loh.calls[loh.calls$Sample.name%in%names(c2)],
                 col = "#1C8FFFFF", data.panel = 1,
                 r0 =at$r0, r1=at$r1, ymax = length(c2))
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 = at$r0, r1=at$r1)
  kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 = at$r0, r1=at$r1)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "LOH", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
  
  # kpRect(kp,data=sample.cn[names(c1)],col = purity.col,chr = "",x0 = -0.05, x1=0, y0 =0,y1 =0.05)
  for(r in seq_len(length(c((length(c2)+7):(total.tracks-6))))){
    rr <- c((length(c2)+7):(total.tracks-6))[r]
    ats <- autotrack(current.track = rr, total.tracks = total.tracks, margin = 0.25)
    sn <- names(c1)[r]
    
    plotCopyNumberCalls(kp, cn.calls = sample.cn[sn] , loh.height = 0.1, loh.column = "loh",  label.cex = 2.5, r0 =ats$r0, r1=ats$r1, cn.colors = cn.colors,  label.margin =0.05)
    
    
    # ats <- autotrack(current.track =rr, total.tracks = total.tracks,margin = 0.15)
    # kpRect(kp,chr="chr1",x0=-150000000,x1 = -100000000, y0 =ats$r0, y1=ats$r1,clipping = FALSE,col =col.purity[names(c1)][r])
    kpRect(kp,chr="chr1",x0=-150000000, x1 = -100000000,r0 = ats$r0, r1 = ats$r1, y0 =0, y1=1, clipping = TRUE, col =col.purity[names(c1)][r], track.margin = 0.25)
    
  }

  at <- autotrack(current.track = c(((total.tracks-6)+1):((total.tracks-6)+2)),total.tracks = total.tracks, margin = 0.05)
  plotCopyNumberSummary(kp, cn.calls = sample.cn[names(c1)] , label.cex = 2.5, labels = "",r0 =at$r1, r1=at$r0,
                        normal.color = "white", gain.color = NA, loss.color = "darkgreen",  data.panel = 1,add.grid = F)
  
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Losses", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 =at$r1, r1=at$r0)
  kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
  
  # kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
  
  at <- autotrack(current.track = c(((total.tracks-6)+3):((total.tracks-6)+4)),total.tracks = total.tracks, margin = 0.65)
  plotCopyNumberSummary(kp, cn.calls = sample.cn[names(c1)] , label.cex = 2.5,labels = "",r0 =at$r1, r1=at$r0,
                        normal.color = "white",loss.color = NA, gain.color ="#FFD700", data.panel = 1)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Gains", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = c("100%", "50%","0"), cex = cex.axis,r0 =at$r1, r1=at$r0)
  kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r1, r1=at$r0)
  
  at <- autotrack(current.track = c(((total.tracks-6)+5):((total.tracks-6)+6)),total.tracks = total.tracks, margin = 0.50)
  kpPlotCoverage(kp, data = loh.calls[loh.calls$Sample.name%in%names(c1)],
                 col = "#1C8FFFFF", data.panel = 1,
                 r0 =at$r0, r1=at$r1, ymax = length(c1))
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 = at$r0, r1=at$r1)
  kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 = at$r0, r1=at$r1)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "LOH", side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 = at$r0, r1=at$r1)
  
  
  dev.off()
}



####################  Biology: explicacion de la tetraploidizacion
groups <- list(C1 = c1,C2= c2)
i=1
cn.prop <- list()
for(i in seq_len(length(groups))){
  clust <- names(groups)[i]
  sn.groups <- groups[[i]]
  sample.cnv.gr <- sample.cn[names(sample.cn) %in% names(sn.groups)]
  c1.cn <- unlist(as(sample.cnv.gr,"GRangesList"))
  c1.cn$width.1 <- width(c1.cn)
  c1.cn <- c1.cn[!seqnames(c1.cn) %in% c("chrX","chrY")]
  c1.cn <- filterChromosomes(c1.cn)
  width.tot.c1 <- sum(c1.cn$width.1) 
  
  for(j in seq_len(length(unique(c1.cn$Sample.name)))){
    sn <- unique(c1.cn$Sample.name)[j]
    c1.sn <-  c1.cn[c1.cn$Sample.name%in%sn]
    
    c1.cn.even <- c1.sn[lapply(c1.sn$cn, "%%", 2) == 0]
    c1.cn.even.width <- sum(c1.cn.even$width.1)
    prop.even.c1 <- c1.cn.even.width/3e9
    
    c1.cn.odd <- c1.sn[!lapply(c1.sn$cn, "%%", 2) == 0]
    c1.cn.odd.width <- sum(c1.cn.odd$width.1)
    prop.odd.c1 <- c1.cn.odd.width/3e9
    
    cn.prop[[sn]] <- data.frame(sample.name = sn, group = clust, even = prop.even.c1*100, odd = prop.odd.c1*100)
  }
  
}
cn.prop <- do.call(rbind,cn.prop)
names(cn.prop)
df.cn.prop <- data.frame(sample.name= rep(cn.prop$sample.name,each =2), groups = rep(cn.prop$group, each=2))
df.cn.prop$EvenOdd <- rep(c("even","odd"),length(df.cn.prop))
df.cn.prop$cn.prop <- as.vector(rbind(cn.prop$even, cn.prop$odd))

svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_CNVprop_cluster1.2_v5.svg")), width = 12, height = 10)

par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")

tt <-boxplot(cn.prop~EvenOdd+groups, data = df.cn.prop,
        col = c( "#7BA7FFFF","#7BA7FFFF","#FFEC6DFF","#FFEC6DFF"),
        frame=FALSE, 
        lwd=3,
        xaxt="n",
        las =2,
        cex.axis =2, yaxt="n",
        ylim = c(0,100),  
        ylab="% Whole genome CN-LOH" ,
        cex.lab = 5, 
        xlab="")
# axis(1, at = c(0:6),labels = rep("",7),lwd.ticks = 0, lwd = 3)

axis(2, at = c(0,50,100), las=2, lwd = 3, cex.axis=5, ylim =c(0,100), ylab="% Odd/Even global CN")

dev.off()
tt <- boxplot(cn.prop~EvenOdd+groups, data = df.cn.prop,
              col = c( "#7BA7FFFF","#7BA7FFFF","#FFEC6DFF","#FFEC6DFF"),
              frame=FALSE,xaxt="n",cex.axis =2, ylim = c(0,100), 
              ylab="" ,cex.lab = 2, las =2 )
dev.off()


# cn.prop.mat <- data.matrix(df.cn.prop)
# cn.prop[cn.prop$EvenOdd =="even",] 
# svg(filename=file.path(execution.dir,"results","images", paste0("Barplot_LOH.DNAprop_cluster1.2_v4.svg")), width = 10, height = 8)
# barplot(df.cn.prop$cn.prop  ,beside = T,ylim = c(0,100), names.arg = c(names(c1),names(c2)), las =2,
#         cex.axis =2, cex.names = 2 )
# dev.off()


#statistical test 
lmd <- lm(cn.prop~groups+EvenOdd+groups:EvenOdd, data= df.cn.prop)
summary(lmd)
tt <- aov(lmd)
TukeyHSD(x =tt,which =  "groups:EvenOdd")
hist(df.cn.prop$cn.prop)

########### loh proportion DNA

loh.prop <- list()
i=1
for(i in seq_len(length(groups))){
  clust <- names(groups)[i]
  sn.groups <- groups[[i]]
  sample.cnv.gr <- sample.cn[names(sample.cn) %in% names(sn.groups)]
  c1.cn <- unlist(as(sample.cnv.gr,"GRangesList"))
  gr.sn <-unique(c1.cn$Sample.name)
  c1.cn$width.1 <- width(c1.cn)
  c1.cn <- c1.cn[!seqnames(c1.cn) %in% c("chrX","chrY")]
  c1.cn <- filterChromosomes(c1.cn)
  width.tot.c1 <- sum(c1.cn$width.1) 
  c1.cn.loh <- c1.cn[c1.cn$loh == TRUE] #For global LOH
  c1.cn.loh <- c1.cn.loh[c1.cn.loh$cn ==2] #For CN-LOH 2n
  j=1
  names(sample.cnv.gr)
  for(j in seq_len(length(gr.sn))){
    sn <- gr.sn[j]
    c1.loh <-  c1.cn.loh[c1.cn.loh$Sample.name%in%sn]
    width.tot<- sum(c1.loh$width.1) 
    
    # ss <- extendRegions(trim(filt.list.sv[[sn]]),extend.start =extend1, extend.end = extend1)
    
    loh.prop[[sn]] <- data.frame(sample.name = sn, groups = clust, loh.prop = width.tot/3e9)
  }
  
  
}

loh.prop.df <- do.call(rbind,loh.prop)
loh.prop.df$loh.prop <- loh.prop.df$loh.prop*100

# png(filename = file.path(execution.dir,"results","images",paste0("BoxPlot_LOH.DNAprop_groups_v3.png")),width = 1200, height = 1000)
svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_CN-LOH.DNAprop_Clusters_v4.svg")), width = 12, height = 10)
par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")

boxplot(loh.prop~groups, data = loh.prop.df, 
        col = c( "#7BA7FFFF","#FFEC6DFF"),
        frame=FALSE, 
        lwd=3,
        xaxt="n",
        las =2,
        cex.axis =2, yaxt="n",
        ylim = c(0,100),  
        ylab="% Whole genome CN-LOH" ,
        cex.lab = 5, 
        xlab="")
axis(1, at = c(0:6),labels = rep("",7),lwd.ticks = 0, lwd = 3)

axis(2, at = c(0,50,100), las=2, lwd = 3, cex.axis=5, ylim =c(0,100), ylab="% Whole genome CN-LOH")
dev.off()

svg(filename=file.path(execution.dir,"results","images", paste0("Barplot_CN-LOH.DNAprop_cluster1.2_v3.svg")), width = 10, height = 8)
barplot(data.matrix(loh.prop.df$loh.prop),beside = T,col = c(rep("#7BA7FFFF", length(c1)), rep("#FFEC6DFF",length(c2))),ylim = c(0,100),names.arg = c(names(c1),names(c2)), las =2,
        cex.axis =2, cex.names = 2 )
dev.off()

#statistical test 
class(loh.prop.df$groups)
loh.prop.df$groups <- as.factor(loh.prop.df$groups)
ttest <- t.test(loh.prop~groups, data= loh.prop.df)
ttest


######Only cn-loh proportion by selected regions

reg.def <- toGRanges(c("chr1:792501-121605200","chr4:69201-190048600","chr9:670035-39500291","chr10:45701-132557735","chr11:194501-135075600",
                       "chr14:80144882-106883718","chr16:46423001-90100500","chr17:404091-35766416"))

width.reg <- sum(width(reg.def))

cn.loh.prop <- list()
i=1
for(i in seq_len(length(groups))){
  clust <- names(groups)[i]
  sn.groups <- groups[[i]]
  sample.cnv.gr <- sample.cn[names(sample.cn) %in% names(sn.groups)]
  c1.cn <- unlist(as(sample.cnv.gr,"GRangesList"))
  gr.sn <-unique(c1.cn$Sample.name)
  c1.cn$width.1 <- width(c1.cn)
  c1.cn <- c1.cn[!seqnames(c1.cn) %in% c("chrX","chrY")]
  c1.cn <- filterChromosomes(c1.cn)
  width.tot.c1 <- sum(c1.cn$width.1) 
  c1.cn.loh <- c1.cn[c1.cn$loh == TRUE]
  c1.cn.loh <- c1.cn.loh[c1.cn.loh$cn ==2]
  c1.cn.loh.sn <- split(c1.cn.loh, c1.cn.loh$Sample.name)
  
  prop.cn.c1.loh.regions <- lapply(c1.cn.loh.sn,function(x){
    y <- subsetByOverlaps(x,reg.def,type = "any",)
    y <- GenomicRanges::intersect(y,reg.def)
    prop <- sum(width(y))/width.reg
    return(prop)
  })
  df <- data.frame("cn.loh.reg"=do.call(rbind,prop.cn.c1.loh.regions))
  # colnames(df) <- "cn.loh.reg"
  df$groups <- clust
  df$Sample.names <- rownames(df)
  cn.loh.prop[[clust]] <- df
}
cn.loh.prop <- do.call(rbind,cn.loh.prop) 
cn.loh.prop$cn.loh.reg <- cn.loh.prop$cn.loh.reg*100

# png(filename = file.path(execution.dir,"results","images",paste0("BoxPlot_specific_CN-LOH.DNApropRegions_groups_v3.png")),width = 1200, height = 1000)
svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_specific_CN-LOH.DNApropRegions_Clusters_v4.svg")), width = 12, height = 10)
par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")

boxplot(cn.loh.reg~groups, data = cn.loh.prop, 
        col = c( "#7BA7FFFF","#FFEC6DFF"),
        frame=FALSE, 
        lwd=3,
        xaxt="n",
        las =2,
        cex.axis =2, yaxt="n",
        ylim = c(0,100),  
        ylab="% CN-LOH region" ,
        cex.lab = 5, 
        xlab="")
axis(1, at = c(0:6),labels = rep("",7),lwd.ticks = 0, lwd = 3)

axis(2, at = c(0,50,100), las=2, lwd = 3, cex.axis=5, ylim =c(0,100), ylab="% CN-LOH regions")
dev.off()

#statistical test
ttest <- t.test(cn.loh.reg~groups, data= cn.loh.prop)
ttest







length(c1.cn.loh)
c1.cn.loh <- c1.cn[c1.cn$loh == TRUE]
c1.cn.loh <-c1.cn.loh[c1.cn.loh$cn==2]
c1.cn.loh.sn <- split (c1.cn.loh, c1.cn.loh$Sample.name)
width.reg <- sum(width(reg.def))
x=""
# x = c1.cn.loh.sn$`MPNST-10_Miller`
prop.cn.c1.loh.regions <- lapply(c1.cn.loh.sn,function(x){
  y <- subsetByOverlaps(x,reg.def,type = "any",)
  y <- intersect(y,reg.def)
  prop <- sum(width(y))/width.reg
  return(prop)
})

c2.cn.loh <- c2.cn[c2.cn$loh == TRUE]
c2.cn.loh <-c2.cn.loh[c2.cn.loh$cn==2]
c2.cn.loh.sn <- split (c2.cn.loh, c2.cn.loh$Sample.name)
prop.cn.c2.loh.regions <- lapply(c2.cn.loh.sn,function(x){
  y <- subsetByOverlaps(x,reg.def,type = "any",)
  y <- intersect(y,reg.def)
  prop <- sum(width(y))/width.reg
  return(prop)
})

list.loh.sn <- c(prop.cn.c1.loh.regions, prop.cn.c2.loh.regions)




cn.loh.df <- data.frame("cn-loh" = unlist(list.loh.sn))
cn.loh.df$sample <- rownames(cn.loh.df)
cn.loh.df <- cn.loh.df[names(sample.names),]
cn.loh.df <- cn.loh.df[!is.na(cn.loh.df$cn.loh),]
sns <- sample.names[sample.names%in% cn.loh.df$sample]
cn.loh.df$sample <- names(sns)
rownames(cn.loh.df) <- names(sns)
cn.loh.df <- cn.loh.df[names(c(c1,c2)),]
cn.loh.df$cluster <- c(rep("c1",length(c1)),rep("c2",length(c2)))
names(sample.cn)
cn.loh.df$cn.loh <- cn.loh.df$cn.loh*100
svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_specific_CN-LOH.DNApropRegions_cluster1.2_v2.svg")), width = 7, height = 8)
tt <- boxplot(cn.loh~cluster, data = cn.loh.df, col = c( "#7BA7FFFF","#FFEC6DFF"),frame=FALSE,xaxt="n",las =2,cex.axis =2, ylim = c(0,100),  ylab="" ,cex.lab = 2, xlab="")
dev.off()
tt <- t.test(cn.loh~cluster,data = cn.loh.df)

cn.loh.mat <- data.matrix(cn.loh.df)
cn.loh.mat <- cn.loh.mat[,1]
cn.loh.mat <- cn.loh.mat[c(names(c1)[length(c1):1],names(c2)[length(c2):1])] 
mean.cn.loh <- mean(cn.loh.mat)
abh <- sd(cn.loh.mat[1:14])
svg(filename=file.path(execution.dir,"results","images", paste0("BarPlot_specific_CN-LOH.DNApropRegions_cluster1.2_v2.svg")), width = 7, height = 8)
barplot(cn.loh.mat,beside = F, las = 2,ylim = c(0,100), col = c(rep("#7BA7FFFF", length(c1)), rep("#FFEC6DFF",length(c2))))
abline(h =mean.cn.loh)
dev.off()
abline(h =mean.cn.loh-abh)
sample.data$
seqnames(reg)
# c1.cn.loh <- subsetByOverlaps(c1.cn.loh, reg.genes)
cn.loh.mat*100
length(cn.loh.mat)
sample.names
dfs <- cbind(loh.prop,cn.loh.prop =cn.loh.df$cn.loh,p.name = cn.loh.df$sample)
dfs.cell <- cbind(loh.prop,cn.loh.prop = cn.loh.df$cn.loh, p.name =cn.loh.df$sample)
dfss <-rbind(dfs,dfs.cell)             
write.table(dfss, file.path("./CN.LOH.prop.RealSet.csv"),col.names = TRUE,row.names = TRUE,sep = "\t")

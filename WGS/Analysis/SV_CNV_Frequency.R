#########################
# SV and CNV frequency  #
#########################
#Packages needed
library(karyoploteR)
library(CopyNumberPlots)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(regioneR)
library(dplyr)
#Parameters

lumpy.file <- "All_Lumpy_SV_Regions_FilteredCommonSVinSamples_mw6.v7.RData"
data.test <- read.table(file = file.path("Supplementary_Table_S8"),sep="\t",header = T)[,c(2:4)]
centromere.gr <- read.table("./Centromeres_hg38_UCSC.txt", sep = "\t", header = T)[,c(2:4)]
centromere.gr <- toGRanges(centromere.gr,genome = "hg38")

execution.dir <- file.path(".")
# setwd(execution.dir)
load(file.path(execution.dir,"results/", lumpy.file))
real.sv <- lapply(filt.list, function(x) x[["AllRegionsFiltered"]])

sample.data <- read.table(file.path(execution.dir,"Sample.info.csv"), sep = "\t", header = T,comment.char = "")
sample.names <- sample.data$Sample.name

sample.data <- sample.data[-c(1,30:32),]
sample.data <- sample.data[sample.data$intoPaper ==TRUE,]
sample.names <- sample.data$Bioinfo.name
names(sample.names) <- sample.data$Sample.name
sample.names <- sample.names[names(sample.names)!=""]
sample.data <- sample.data[sample.data$Bioinfo.name %in% sample.names,]
pur <- sample.data$Tumor.purity.est

sns <- sample.names
for(i in seq_len(length(sns))){
  sn <- sns[i]
  names(real.sv)[names(real.sv)%in%names(sn)] <-sn
}

execution.dir <- file.path(".")
# setwd(execution.dir)
load(file.path(execution.dir,"results/", lumpy.file))
val.sv <- lapply(filt.list, function(x) x[["AllRegionsFiltered"]])

sample.data <- read.table(file.path(execution.dir,"Sample.info.1.csv"), sep = ",", header = TRUE)
sample.names <- sample.data$Sample.name

sample.data <- sample.data[sample.data$BAM_Tumor ==TRUE,]
sample.names <- sample.data$Sample.name
names(sample.names) <- sample.data$File.name
# sns <- c()
# sns <- c(sns,sample.names,"ST88-14"="ST88-14","STS-26T"="STS-26T","HS-Sch-2"="HS-Sch-2","HS-PSS"="HS-PSS","sNF96.2"="sNF96.2","NF90-8"="NF90-8","S462"="S462","NMS-2"="NMS-2")
sns <- c(sns,sample.names) #no lines
SVs <- c(real.sv, val.sv)

scoring.table <- read.table("Supplementary_Table_S8", sep ="\t",header = T)
#NoMPNST cell lines
scoring.table <- scoring.table[grepl("MPNST",scoring.table$Sample.names),]

SVs <- SVs[names(SVs)%in%scoring.table$Sample.name]
names(SVs)
names(SVs[scoring.table$Sample.name])
# #only Miller Samples
# scoring.table$Groups <-"0" 
# scoring.table$Groups[1:38] <- "1"
# scoring.table$Groups[47:53] <- "2"
# scoring.table$Groups[54:59] <- "3"

scoring.table <- scoring.table[scoring.table$Sample.name %in% names(SVs),]
# scoring.table$Groups_order[scoring.table$Groups_order==5] <-4
SVs <- SVs[names(SVs)%in%scoring.table$Sample.name]

groups <- c(1,2,3,4)
names(groups) <- c("#7C2C00",  "tomato","goldenrod", "#E09B9B")
groups
g1.sn <- scoring.table$Sample.name[scoring.table$Groups == 1]
g2.sn <- scoring.table$Sample.name[scoring.table$Groups == 2]
g3.sn <- scoring.table$Sample.name[scoring.table$Groups == 3]
g0.sn <- scoring.table$Sample.name[scoring.table$Groups == 4]

# Percentage of SV affecting the genome
filt.list.sv <- SVs
filt.list.sv<- lapply(filt.list.sv,function(x){
  names.todel <- names(subsetByOverlaps(x,centromere.gr, type = "within"))
  x[!names(x)%in%names.todel]
})

# filt.list.sv <- lapply(filt.list.sv, function(x) x[["AllRegionsFiltered"]])
# filt.list.sv <- lapply(filt.list.sv, function(x) x+1e6)
red.regions.svs <- list()

for(gr in seq_len(length(groups))){
  g <- groups[gr]
  g.sn <- scoring.table$Sample.name[scoring.table$Groups_order == g]
  filt.list.sv <- SVs
  filt.list.sv <- lapply(filt.list.sv, function(x) x+1e6)
  
  sv.g.sn <- filt.list.sv[names(filt.list.sv) %in% g.sn] 
  sv.g.sn <- lapply(sv.g.sn, function(x)reduce(x))
  sv.g.sn <- unlist(as(sv.g.sn,"GRangesList"))
  
  svs <- c("DUP", "DEL","INV","ICR")
  sv.red <- list()
  for(i in seq_len(length(svs))){
    sv <- svs[i] 
    sv.red[[sv]]  <- lapply(filt.list.sv, function(x) reduce(x[x$SVTYPE==sv],drop.empty.ranges = T))
  } 
  
  
  filt.list.sv <- unlist(as(filt.list.sv,"GRangesList"))
  
  for(j in seq_len(length(svs))){
    sv <- svs[j]
    regions <- filt.list.sv[filt.list.sv$SVTYPE == sv] 
    red.reg <- unlist(as(sv.red[[sv]],"GRangesList"))
    
    red.reg$SharedBy<- NA
    red.reg.gr <- red.reg[names(red.reg)%in% g.sn] 
    for(i in seq_len(length(red.reg.gr))){
      
      reduced.regions <-subsetByOverlaps(regions[regions$Sample %in% g.sn] ,red.reg.gr[i])
      red.reg.gr$SharedBy[i] <- length(unique(reduced.regions$Sample))
    }
    red.regions.svs[[paste0("group_",g,"_",sv)]] <- red.reg.gr
    
  }
  red.regions.svs[[paste0("group_",g,"_","AllSV")]] <- sv.g.sn
  
}


#Porportion of DNA affected by SVs
sns <- sns[sns %in% scoring.table$Sample.name]
# filt.list.icr.inv <- lapply(SVs, function(x) x[x$SVTYPE %in% c("INV","ICR")])
filt.list.icr.inv <- SVs
filt.list.icr.inv<- lapply(filt.list.icr.inv,function(x){
  names.todel <- names(subsetByOverlaps(x,centromere.gr, type = "within"))
  x[!names(x)%in%names.todel]
})
#cantidad de genoma afectada por SVs
extend1 <- 1e6
svs.mb <- list()
for(i in seq_len(length(sns))){
  sn <- sns[i]
  if(!sn %in% names(filt.list.icr.inv))next
  ss <- extendRegions(trim(filt.list.icr.inv[[sn]]),extend.start =extend1, extend.end = extend1)
  svs.mb[[sn]] <- reduce(ss)
}
svs.mb <-svs.mb[names(svs.mb) %in%scoring.table$Sample.name]
res <- unlist(lapply(svs.mb,function(x){sum(width(x))}))/3e9
round(res,digits = 3)
col.gr <- scoring.table$Groups_order
for(i in seq_len(length(groups))){
  gr <- groups[i]
  col.gr[col.gr %in% gr] <- names(gr)
}
barplot(res, las =2,col = col.gr)

#mean by group of SVs
mean.gr.sv <- list()
for(i in seq_len(length(groups))){
  gr  <- groups[i]
  gr.sn <- scoring.table$Sample.name[scoring.table$Groups_order%in%gr]
  mean.gr.sv[[paste0("group_",gr)]] <- mean(res[names(res)%in%gr.sn])
}
mean.gr.sv <- do.call(cbind,mean.gr.sv)
barplot(mean.gr.sv)
res <- res[scoring.table$Sample.names]
data.res.sv <- data.frame(res, groups = factor(scoring.table$Groups_order,levels=c(1,2,3,4)))
# data.res <- data.frame(res, groups = factor(scoring.table$Cluster,levels=c("C1","C2")))
data.res.sv <- data.res.sv[!is.na(data.res.sv$groups),]

png(filename = file.path(execution.dir,"results","images",paste0("Boxplot_Freq_SVs_GenomicGroupSamples_V6.png")),width = 1200, height = 1000)
svg(filename = file.path(execution.dir,"results","images",paste0("Boxplot_Freq_SVs_GenomicGroupSamples_V7.svg")),width = 12, height = 10)

par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")
bp <- boxplot(data.res.sv$res~data.res.sv$groups,
              lwd=3,
              xlab = "",
              ylab="% SVs",
              main = "",
              boxfill=names(groups), 
              frame=FALSE,
              cex.axis = 5,
              cex.main =3,
              cex.lab = 5,
              xaxt="n",yaxt="n", ylim= c(0, 1 ))

axis(1, at = c(0:6),labels = rep("",7),lwd.ticks = 0, lwd = 3)


axis(2, at = c(0,0.5,1),labels = c("0","50","100"), las=2, lwd = 3, cex.axis=5, ylim =c(0,1))

# legend("topright", inset=c(-0.4,0),box.lwd = "", xpd = TRUE, fill = unique(names(groups)), legend = groups, horiz = F, cex = 3,col = 2)

dev.off()


boxplot(data.res.sv$res~data.res.sv$groups,col= names(groups) )

lmd <- lm(data.res.sv$res~data.res.sv$groups)
summary(lmd)
tt <- aov(lmd)
TukeyHSD(x =tt)
# scoring.table$SV.freq <- c()
# scoring.table$SV.req[which(scoring.table$Sample.name%in%names(res))]$SV.req <- res
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
cn.colors <- transparent(cn.colors)
chrs <- c("canonical",as.character(seqlevels(filt.list.sv)))

chrs <- "canonical"
label.axis <- c("0","50%","100%")
ymax = 1
cex.axis = 2
cex.titel.axis=3
svs <- c(svs,"AllSV")
gr <- 4
for(gr in seq_len(length(groups))){
  g <- groups[gr]
  g.sn <- scoring.table$Sample.name[scoring.table$Groups == g]
  col.gr <- names(g)
  
  for(k in seq_len(length(chrs))){
    chr <- chrs[k] 
    
    svg(filename=file.path(execution.dir,"results","images", paste0(chr,"_SV_group_",g,"_Clasifier_Reduced_V12.svg")), width = 30, height = 20)
    # png(filename=file.path(execution.dir,"results","images", paste0(chr,"_SV_cluster1.2_ReducedV4.png")), width = 2000, height = 1800)
    
    kp <- karyoploteR::plotKaryotype(genome = "hg38", plot.type = 4, chromosomes = chr, plot.params = pp,ideogram.plotter = NULL, labels.plotter = NULL)
    karyoploteR::kpAddCytobandsAsLine(kp)
    kpAddChromosomeNames(kp, cex = 2.5, srt = 45 )
    # kpPlotMarkers(kp, data = tsgs, labels = tsgs$SYMBOL, cex =1, lwd = 1, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.05, 0.002),ignore.chromosome.ends = TRUE,clipping = FALSE)
    j = 2
    for(j in seq_len(length(svs))){
      sv <- svs[j]
      
      
      at <- autotrack(current.track = j, total.tracks = 5,margin = 0.25)
      
      # kpRect(kp, data = red.regions.svs[[paste0("c1.",sv)]] ,x0 = start(red.regions.svs[[paste0("c1.",sv)]] ), x1 = end(red.regions.svs[[paste0("c1.",sv)]] ),y0 =0,y1 =red.regions.svs[[paste0("c1.",sv)]] $SharedBy,ymax = length(c1),
      #        col = transparent("black"), data.panel = 1,
      #        r0 =at$r0, r1=at$r1)
      kpPlotCoverage(kp,data=red.regions.svs[[paste0("group_",g,"_",sv)]],ymax =length(g.sn),col = col.gr, data.panel = 1,r0 =at$r0, r1=at$r1,border = "grey")
      names(red.regions.svs)
      kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis,r0 =at$r0, r1=at$r1)
      kpAbline(kp, h =c(0,0.5,1), ymax =ymax ,r0 =at$r0, r1=at$r1)
      karyoploteR::kpAddLabels(karyoplot = kp,labels = sv, side = "left", cex=cex.titel.axis,  data.panel = 1,label.margin = 0.05,r0 =at$r0, r1=at$r1)
      
    }
  }
  
  dev.off()
}



################# proportion of DNA affected by SV
filt.list.sv <- SVs
filt.list.sv<- lapply(filt.list.sv,function(x){
  names.todel <- names(subsetByOverlaps(x,centromere.gr, type = "within"))
  x[!names(x)%in%names.todel]
})
# subsetByOverlaps(filt.list.sv$JSL_MPNST,gene.markers.gf[gene.markers.gf$Genes == "CDKN2A"])
#All SVs per group
svs.groups <-list()
for(gr in seq_len(length(groups))){
  g <- groups[gr]
  g.sn <- scoring.table$Sample.name[scoring.table$Groups == g]
  col.gr <- names(g)
  svs.gr <- filt.list.sv[names(filt.list.sv)%in% g.sn]
  svs.gr <- unlist(as(svs.gr,"GRangesList"))
  svs.gr <- svs.gr+1e6
  
  sv.per.sn <- c(split(svs.gr, svs.gr$Sample))
  sv.per.sn <- unlist(lapply(sv.per.sn, length))
  svs.gr <- split(svs.gr, svs.gr$SVTYPE)
  sv.type <- lapply(svs.gr,function(x)table(x$Sample))
  sv.type <- lapply(sv.type,function(x) x/(sv.per.sn[names(x)]) )
  
  sv.type <- lapply(sv.type,function(x)sum(x)/length(g.sn))
  svs.gr <- unlist(sv.type)
  svs.groups[[paste0("SVs_group",g)]] <- svs.gr
}
count  <- do.call(cbind, svs.groups)
for(i in seq_len(ncol(count))){
  col.count <- colnames(count)[i]
  count[,col.count] <- (count[,col.count]/ sum(count[,col.count]))*100
  
}


svg(filename=file.path(execution.dir,"results","images", paste0("All.SV.1mb.Groups.Sample.prop_V7.svg")), width = 5, height = 10)
barplot(count,
        col = c("#4ED058","#FFA300","grey80","darkblue"),
        legend.text = c("DEL","DUP", "Inter-CR","Intra-CR"),ylim=c(0,100),main= "All SVs",
        las=2,cex.axis = 2,cex.names=2.5)

dev.off()


############ CNV proportions #########################
# scoring.table <- read.table("./ContingencyTable/ScoringAndInformationTable.NewScoringORDER.wholeinfo.v12.csv", sep ="\t",header = T)
load(file.path("./CN_AllMPNSTpaper_LOHcorrected_RealValidationSets.RData"))
groups <- c(1,2,3,4)
names(groups) <- c("#7C2C00",  "tomato","goldenrod", "#E09B9B")
# save(sample.cn, file = file.path("./CN_AllMPNSTpaper_LOHcorrected_RealValidationSets.RData"))

# sn.cn$Sample.name <- sn
# mcols(sn.cn) <- mcols(sn.cn)[c(1:10,length(mcols(sn.cn)),11,12)]
# sample.cn[["MPNST_54_B"]] <- sn.cn
# names(sample.cn)
sample.cn <- sample.cn[names(sample.cn)%in%scoring.table$Sample.name]
names(sample.cn)
# names(sample.cn)[!names(sample.cn)%in%scoring.table$Sample.name]
sample.cn<- lapply(sample.cn,function(x){
  names.todel <- names(subsetByOverlaps(x,centromere.gr, type = "within"))
  x[!names(x)%in%names.todel]
})
cn.prop <- sample.cn
for(i in seq_len(length(cn.prop))){
  sn <- names(cn.prop)[i]
  x <- cn.prop[[sn]]
  x$cn[x$cn==1]<- 1
  x$cn[x$cn<1]<- 0
  x$cn[x$cn==2]<- 2
  x$cn[x$cn>2]<- 3
  cn.prop[[sn]] <-x
}


#cantidad de genoma afectada por CNVs
extend1 <- 1e6
names(sns)<-sns
sns <- sns[scoring.table$Sample.name]
svs.mb <- list()
for(i in seq_len(length(sns))){
  sn <- sns[i]
  x <- cn.prop[[sn]]
  j=1
  for(j in seq_len(length(c(0,1,2,3)))){
    cnv <- c(0,1,2,3)[j]
    cn <-  x[x$cn==cnv]
    ss <- extendRegions(trim(cn),extend.start =extend1, extend.end = extend1)
    svs.mb[[paste0(sn,"_CNV_freq_",cnv)]] <- reduce(ss)
  }
  
}

res <- unlist(lapply(svs.mb,function(x){sum(width(x))}))/3e9
round(res,digits = 3)
col.gr <- scoring.table$Groups_order
for(i in seq_len(length(groups))){
  gr <- groups[i]
  col.gr[col.gr %in% gr] <- names(gr)
}
barplot(res[grepl("freq_3",names(res))]*100, las =2,col = col.gr)

data.res <- data.frame(res, groups = rep(factor(scoring.table$Groups_order,levels = c(1,2,3,4)),each =4), cnv = factor(rep(c(0,1,2,3)), levels=c(0,1,2,3)))
# data.res <- data.frame(res, groups = factor(scoring.table$Cluster,levels=c("C1","C2")),cnv = factor(rep(c(0,1,2,3)), levels=c(0,1,2,3)))
data.res <- data.res[!is.na(data.res$groups),]

cn.df <- data.res[data.res$cnv ==3,]


png(filename = file.path(execution.dir,"results","images",paste0("Boxplot_Freq_CNV_Gain_GenomicGroupSamples_v4.png")),width = 1200, height = 1000)
svg(filename = file.path(execution.dir,"results","images",paste0("Boxplot_Freq_CNV_Gain_GenomicGroupSamples_v4.svg")),width = 12, height = 10)

  par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")
  bp <- boxplot(cn.df$res~cn.df$groups,
                lwd=3,
                xlab = "",
                ylab="Freq CNV Gain",
                main = "",
                boxfill=names(groups), 
                frame=FALSE,
                cex.axis = 5,
                cex.main =3,
                cex.lab = 5,
                xaxt="n",yaxt="n", ylim= c(0, 1 ))
  
  axis(1, at = c(0:6),labels = rep("",7),lwd.ticks = 0, lwd = 3)

  
  axis(2, at = c(0,0.5,1), las=2, lwd = 3, cex.axis=5, ylim =c(0,1))
  
  legend("topright", inset=c(-0.4,0),box.lwd = "", xpd = TRUE, fill = unique(names(groups)), legend = groups, horiz = F, cex = 3,col = 2)

  dev.off()




boxplot(cn.df$res~cn.df$groups,col= names(groups) )
lmd <- lm(cn.df$res~cn.df$groups)
summary(lmd)
tt <- aov(formula = cn.df$res~cn.df$groups)
TukeyHSD(x =tt)

cnv.freq.0 <-data.res[data.res$cnv ==0,]
cnv.freq.1 <- data.res[data.res$cnv ==1,]
cnv.freq.2 <- data.res[data.res$cnv ==2,]
cnv.freq.3 <- data.res[data.res$cnv ==3,]

# df.scoring <- data.frame(Sample.name =rownames(data.res.sv), SV.freq = data.res.sv$res, Group = data.res.sv$groups, CNV.freq.2 = cnv.freq.2$res, CNV.freq.3 = cnv.freq.3$res)
# scoring.table <- read.table("./ContingencyTable/ScoringAndInformationTable.NewScoringORDER.wholeinfo.v12.csv", sep ="\t",header = T)
# scoring.table$SV.freq[which(scoring.table$Sample.name%in%rownames(data.res.sv))] <- data.res.sv$res
# scoring.table$CNV.freq.2[which(scoring.table$Sample.name%in%df.scoring$Sample.name)] <- cnv.freq.2$res
# scoring.table$CNV.freq.3[which(scoring.table$Sample.name%in%df.scoring$Sample.name)] <- cnv.freq.3$res
# nrow(scoring.table)
# is.null(scoring.table$CNV.freq.2)
# scoring.table$CNV.freq.2[unlist(lapply(scoring.table$CNV.freq.2, function(x) is.null(x)))] <- NA
# scoring.table$CNV.freq.2 <-unlist(scoring.table$CNV.freq.2)
# colnames(scoring.table)
# write.table(x = scoring.table,file = file.path("./ContingencyTable/ScoringAndInformationTable.NewScoringORDER.wholeinfo.v11.csv"),sep = "\t", col.names = TRUE, row.names = FALSE)


n.snvs <- scoring.table$Number.SNVs
names(n.snvs) <- scoring.table$Groups
n.snvs <- data.frame(n.snvs, groups =factor(scoring.table$Groups_order,levels = c( 1,2,3,4) ))
col.gr <- as.character(n.snvs$groups)
for(i in seq_len(length(groups))){
  gr <- groups[i]
  col.gr[col.gr %in% gr] <- names(gr)
}
png(filename = file.path(execution.dir,"results","images",paste0("Barplot_num.SNVs_GenomicGroupSamples_V2.png")),width = 1200, height = 1000)
par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")
barplot(data.matrix(n.snvs$n.snvs),beside = T,col= col.gr,ylab = "# SNVs",las = 2, cex.axis = 5, cex.lab =5)
legend("topright", inset=c(-0.4,0),box.lwd = "", xpd = TRUE, fill = unique(names(groups)), legend = groups, horiz = F, cex = 3,col = 2)

dev.off()

boxplot(n.snvs$n.snvs~n.snvs$groups,col= names(groups) )
lmd <- lm(n.snvs$n.snvs~n.snvs$groups)
summary(lmd)
tt <- aov(lmd)
TukeyHSD(x =tt)



####################  Biology: explanation of tetraploidization################
sn.groups <- scoring.table$Sample.names
names(sn.groups)<-scoring.table$Groups_order
names(sn.groups)[names(sn.groups)==5] <-4

cn.prop <- list()
for(i in seq_len(length(groups))){
  gr <- groups[i]
  sample.cnv.gr <- sample.cn[names(sample.cn) %in% sn.groups[names(sn.groups)==gr]]
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

    cn.prop[[sn]] <- data.frame(sample.name = sn, group = gr, even = prop.even.c1*100, odd = prop.odd.c1*100)
  }
  
}
cn.prop <- do.call(rbind,cn.prop)
df.cn.prop <- data.frame(sample.name= rep(cn.prop$sample.name,each =2), groups = rep(cn.prop$group, each=2))
df.cn.prop$EvenOdd <- rep(c("even","odd"),length(df.cn.prop))
df.cn.prop$cn.prop <- as.vector(rbind(cn.prop$even, cn.prop$odd))
df.cn.prop$color <-""
for(i in seq_len(length(groups))){
  df.cn.prop$color[df.cn.prop$groups == groups[i]] <-names(groups)[i]
}


svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_CNVprop_GenomicGroups_v3.svg")), width = 7, height = 8)
tt <- boxplot(cn.prop~EvenOdd+groups, data = df.cn.prop,
              col = rep(names(groups),each=2),
              frame=FALSE,xaxt="n",cex.axis =2, ylim = c(0,100), 
              ylab="" ,cex.lab = 2, las =2 )
dev.off()

#statistical test 
lmd <- lm(cn.prop~groups+EvenOdd, data= df.cn.prop)
summary(lmd)
tt <- aov(lmd)
tt
TukeyHSD(x =tt,which =  "EvenOddood")



cn.prop.mat <- data.matrix(cn.prop)
cn.prop[cn.prop$EvenOdd =="even",] 
svg(filename=file.path(execution.dir,"results","images", paste0("Barplot_LOH.DNAprop_cluster1.2_v2.svg")), width = 10, height = 8)
barplot(cn.prop$cn.prop  ,beside = T,ylim = c(0,1),names.arg = c(names(c1),names(c2)), las =2,
        cex.axis =2, cex.names = 2 )
dev.off()



########### loh proportion DNA ####################
loh.prop <- list()
i=1
for(i in seq_len(length(groups))){
  gr <- groups[i]
  sample.cnv.gr <- sample.cn[names(sample.cn) %in% sn.groups[names(sn.groups)==gr]]
  c1.cn <- unlist(as(sample.cnv.gr,"GRangesList"))
  gr.sn <-unique(c1.cn$Sample.name)
  c1.cn$width.1 <- width(c1.cn)
  c1.cn <- c1.cn[!seqnames(c1.cn) %in% c("chrX","chrY")]
  c1.cn <- filterChromosomes(c1.cn)
  width.tot.c1 <- sum(c1.cn$width.1) 
  c1.cn.loh <- c1.cn[c1.cn$loh == TRUE]
  # c1.cn.loh <- c1.cn.loh[c1.cn.loh$cn ==2]
  j=1
  names(sample.cnv.gr)
  for(j in seq_len(length(gr.sn))){
    sn <- gr.sn[j]
    c1.loh <-  c1.cn.loh[c1.cn.loh$Sample.name%in%sn]
    width.tot<- sum(c1.loh$width.1) 
    
    # ss <- extendRegions(trim(filt.list.sv[[sn]]),extend.start =extend1, extend.end = extend1)
    
    loh.prop[[sn]] <- data.frame(sample.name = sn, groups = gr, loh.prop = width.tot/3e9)
  }
  
  
}

loh.prop.df <- do.call(rbind,loh.prop)
loh.prop.df$loh.prop <- loh.prop.df$loh.prop*100

png(filename = file.path(execution.dir,"results","images",paste0("BoxPlot_LOH.DNAprop_groups_v4.png")),width = 1200, height = 1000)
svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_LOH.DNAprop_groups_v6.svg")), width = 12, height = 10)
par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")

boxplot(loh.prop~groups, data = loh.prop.df, 
              col = names(groups),
              frame=FALSE, 
              lwd=3,
              xaxt="n",
              las =2,
              cex.axis =2, yaxt="n",
              ylim = c(0,100),  
              ylab="% Whole genome LOH" ,
              cex.lab = 5, 
              xlab="")
axis(1, at = c(0:6),labels = rep("",7),lwd.ticks = 0, lwd = 3)

axis(2, at = c(0,50,100), las=2, lwd = 3, cex.axis=5, ylim =c(0,100), ylab="% Whole genome CN-LOH")
dev.off()
#statistical test 
class(loh.prop.df$groups)
loh.prop.df$groups <- as.factor(loh.prop.df$groups)
lmd <- lm(loh.prop~groups, data= loh.prop.df)
summary(lmd)
tt <- aov(lmd)
TukeyHSD(x =tt,which =  "groups")

######Only cn-loh proportion
reg.def <- toGRanges(c("chr1:792501-121605200","chr4:69201-190048600","chr9:670035-39500291","chr10:45701-132557735","chr11:194501-135075600",
                       "chr14:80144882-106883718","chr16:46423001-90100500","chr17:404091-35766416"))

width.reg <- sum(width(reg.def))

cn.loh.prop <- list()
i=1
for(i in seq_len(length(groups))){
  gr <- groups[i]
  sample.cnv.gr <- sample.cn[names(sample.cn) %in% sn.groups[names(sn.groups)==gr]]
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
  df$groups <- gr
  df$Sample.names <- rownames(df)
  cn.loh.prop[[gr]] <- df
}
cn.loh.prop <- do.call(rbind,cn.loh.prop) 
cn.loh.prop$cn.loh.reg <- cn.loh.prop$cn.loh.reg*100

png(filename = file.path(execution.dir,"results","images",paste0("BoxPlot_specific_CN-LOH.DNApropRegions_groups_v3.png")),width = 1200, height = 1000)
svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_specific_CN-LOH.DNApropRegions_groups_v.svg")), width = 12, height = 10)
par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")

boxplot(cn.loh.reg~groups, data = cn.loh.prop, 
        col = names(groups),
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

axis(2, at = c(0,50,100), las=2, lwd = 3, cex.axis=5, ylim =c(0,100), ylab="% Whole genome CN-LOH")
dev.off()


  
svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_specific_CN-LOH.DNApropRegions_groups_v3.svg")), width = 7, height = 8)
tt <- boxplot(cn.loh.reg~groups, data = cn.loh.prop, col = names(groups),frame=FALSE,xaxt="n",las =2,cex.axis =2, ylim = c(0,100),  ylab="% CN-LOH regions" ,cex.lab = 2, xlab="")
dev.off()
  
#statistical test 
lmd <- lm(cn.loh.reg~groups, data= cn.loh.prop)
summary(lmd)
tt <- aov(lmd)
TukeyHSD(x =tt,which =  "groups")



######## Chr1 Unbalance: 1p LOH and CN<1q ####
kp <- plotKaryotype(genome = "hg38",chromosomes = "chr1")
plotCopyNumberCalls(kp,cn.calls = x)
dd <- karyoploteR::getCytobands(genome = "hg38")
dd <- filterChromosomes(dd)
dd <- dd[seqnames(dd)=="chr1"]
chr1p <- reduce(dd[grepl("p",dd$name)])
chr1q <- reduce(dd[grepl("q",dd$name)])
# x= sample.cn.groups[[3]]
sample.cn.groups <- sample.cn[names(sample.cn)%in% sn.groups]
names(sample.cn.groups)
chr1.list <- lapply(sample.cn.groups, function(x){
  p <- subsetByOverlaps(x,chr1p)
  if(length(p)==1){
   p <- p[grepl("p", q$chr.band)]
    
  }else{
    p <- p[!grepl("q", p$chr.band)]
    
  }
  # p <- p[!grepl("q", p$chr.band)]
  
  p.cn <- data.frame(cn.p =p$cn)
  p.cn$width <- width(p)
  width.cn <- tapply(p.cn$width, p.cn$cn, FUN=sum)
  mainp.cn <- names(width.cn)[which.max(width.cn)]
  
  p.loh <- data.frame(loh =p$loh)
  p.loh$width <- width(p)
  width.loh <- tapply(p.loh$width, p.loh$loh, FUN=sum)
  mainp.loh <- names(width.loh)[which.max(width.loh)]
  
  
  q <- subsetByOverlaps(x,chr1q)
  if(length(q)==1){
    q <- q[grepl("q", q$chr.band)]
    
  }else{
    q <- q[!grepl("p", q$chr.band)]
    
  }
  q.cn <- data.frame(cn.q =q$cn)
  q.cn$width <- width(q)
  width.cn <- tapply(q.cn$width, q.cn$cn, FUN=sum)
  mainq.cn <- names(width.cn)[which.max(width.cn)]
  
  q.loh <- data.frame(loh =q$loh)
  q.loh$width <- width(q)
  width.loh <- tapply(q.loh$width, q.loh$loh, FUN=sum)
  mainq.loh <- names(width.loh)[which.max(width.loh)]
  df <- data.frame(cn.p= mainp.cn, loh.p = mainp.loh,cn.q= mainq.cn, loh.q = mainq.loh)
})
class(chr1p.list)
df.chr1 <- do.call(rbind,chr1.list)
df.chr1$Groups <- names(sn.groups)
df.chr1$Unbalance <- FALSE
df.chr1$Unbalance[df.chr1$cn.p>1 & df.chr1$loh.p == TRUE & df.chr1$cn.q>df.chr1$cn.p] <- TRUE
# df.chr1$Unbalance[df.chr1$cn.p == 2 & df.chr1$loh.p == TRUE & df.chr1$cn.q>2] <- TRUE
table(df.chr1$Unbalance)
rownames(df.chr1)[df.chr1$Unbalance][21]

df.chr1[rownames(df.chr1)[df.chr1$Unbalance][21],]
scoring.table$Sample.names==rownames(df.chr1)
scoring.table$chr1.p.cn[df.chr1$cn.p == 2 & df.chr1$loh.p == FALSE] <- "2n"
scoring.table$chr1.p.cn[df.chr1$cn.p >2 ] <- "gain"
scoring.table$chr1.p.cn[df.chr1$cn.p ==1] <- "het.loss"
scoring.table$chr1.p.cn[df.chr1$cn.p ==0] <- "hom.loss"
scoring.table$chr1.p.cn[df.chr1$cn.p >1 & df.chr1$loh.p == TRUE] <- "cn.loh"
scoring.table$chr1.q.cn[df.chr1$cn.q >1 & df.chr1$loh.q == TRUE&df.chr1$cn.q>df.chr1$cn.p] <- "cn.loh"

scoring.table$chr1.q.cn[df.chr1$cn.q>2] <- "gain"
scoring.table$chr1.q.cn[df.chr1$cn.q == 2 & df.chr1$loh.q == FALSE] <- "2n"
scoring.table$chr1.q.cn[df.chr1$cn.q == 1] <- "het.loss"
scoring.table$chr1.q.cn[df.chr1$cn.q == 0] <- "hom.loss"
scoring.table$chr1.q.cn[df.chr1$cn.q >1 & df.chr1$loh.q == TRUE] <- "cn.loh"
scoring.table$chr1.q.cn[df.chr1$cn.q >1 & df.chr1$loh.q == TRUE&df.chr1$cn.q>df.chr1$cn.p] <- "cn.loh"

scoring.table$chr1.unbalance <- FALSE
scoring.table$chr1.unbalance[scoring.table$Sample.names %in% rownames(df.chr1)[df.chr1$Unbalance]] <- TRUE

table(scoring.table$chr1.unbalance)













########## Frequency of chrs 2 7 8q and 12 gain
i=1
gr <- groups[i]

sample.cnv.gr <- sample.cn[names(sample.cn) %in% sn.groups[names(sn.groups)==gr]]
chrs.gain <- svs.mb[grepl("CNV_freq_3", names(svs.mb))]

chrs.gain <- lapply(sample.cnv.gr,function(x) x[seqnames(x) %in% c("chr2","chr7","chr8","chr12")])
chrs.gain <- lapply(chrs.gain,function(mcols(x)) mcols(x)$cn <-3)
x = chrs.gain$MPNST_51_CNV_freq_3
mcols(x)$cn <- 3

kp <-plotKaryotype(genome = "hg38", chromosomes = c("chr2","chr7","chr8","chr12"), plot.type = 3)
at <- autotrack(current.track = c(length(g3)+2),total.tracks = total.tracks, margin = 0.65)
plotCopyNumberSummary(kp, cn.calls = sample.cn[rownames(g3)] , label.cex = 2.5,labels = "",r0 =1, r1=0,
                      normal.color = "white",loss.color = NA, gain.color ="#FFD700",data.panel = 1)
karyoploteR::kpAddLabels(karyoplot = kp,labels = "Gains", side = "left", cex=2,  data.panel = 1,label.margin = 0.05,r0 = 0, r1=1)

pp <-plotCopyNumberSummary(karyoplot = kp, cn.calls = sample.cn[rownames(g3)],)
pp$
chrs.gain$MPNST_31_CNV_freq_3[seqnames(chrs.gain$MPNST_31_CNV_freq_3) %in% c("chr2","chr7","chr8","chr12")]
g3

kp <-plotKaryotype(genome = "hg38", chromosomes = c("chr8","chr12"), plot.type = 3,)
plotCopyNumberCalls(kp,cn.calls = sample.cn[rownames(g3)],data.panel =2)
# at <- autotrack(current.track = c(length(g3)+2),total.tracks = total.tracks, margin = 0.65)
plotCopyNumberSummary(kp, cn.calls = sample.cn[rownames(g3)] , label.cex = 2.5,labels = "",r0 =1, r1=0,
                      normal.color = "white",loss.color = NA, gain.color ="#FFD700",data.panel =1 )





coverage(sample.cn[rownames(g3)][mcols(sample.cn$MPNST_NF1_13)])



chrs = c("chr2","chr7","chr8","chr12")
chrs.df <- data.frame(row.names = sn.groups)

i=3
for(i in seq_len(length(chrs))){
  chr <- chrs[i]
  # kp <- plotKaryotype(genome = "hg38",chromosomes = "chr2")
  # plotCopyNumberCalls(kp,cn.calls = x)
  dd <- karyoploteR::getCytobands(genome = "hg38")
  dd <- filterChromosomes(dd)
  dd <- dd[seqnames(dd)==chr]
  chrp <- reduce(dd[grepl("p",dd$name)])
  chrq <- reduce(dd[grepl("q",dd$name)])
  # p <- subsetByOverlaps(sample.cn.groups$MPNST_NF1_13, chrp)
  
  sample.cn.groups <- sample.cn[names(sample.cn)%in% sn.groups]
  names(sample.cn.groups)
  # s.cn= sample.cn.groups$`MPNST-50_Miller`
  chr.list <- lapply(sample.cn.groups, function(s.cn){
    p <- subsetByOverlaps(s.cn,chrp)
    if(length(p)==1){
      p <- p[grepl("p", p$chr.band)]
      
    }else{
      p <- p[!grepl("q", p$chr.band)]
      
    }
    # p <- p[!grepl("q", p$chr.band)]
    if(length(p) == 0){
      mainp.cn <- NA
      mainp.loh <- NA
    }else{
      p.cn <- data.frame(cn.p =p$cn)
      p.cn$width <- width(p)
      width.cn <- tapply(p.cn$width, p.cn$cn, FUN=sum)
      mainp.cn <- names(width.cn)[which.max(width.cn)]
      
      # p.loh <- data.frame(loh =p$loh)
      # p.loh$width <- width(p)
      # width.loh <- tapply(p.loh$width, p.loh$loh, FUN=sum)
      # mainp.loh <- names(width.loh)[which.max(width.loh)]
    }
    
    
    
    q <- subsetByOverlaps(s.cn,chrq)
    if(length(q)==1){
      q <- q[grepl("q", q$chr.band)]
      
    }else{
      q <- q[!grepl("p", q$chr.band)]
      
    }
    
    if(length(q) == 0){
      mainq.cn <- NA
      mainq.loh <- NA
    }else{
      q.cn <- data.frame(cn.q =q$cn)
      q.cn$width <- width(q)
      width.cn <- tapply(q.cn$width, q.cn$cn, FUN=sum)
      mainq.cn <- names(width.cn)[which.max(width.cn)]
      
      # q.loh <- data.frame(loh =q$loh)
      # q.loh$width <- width(q)
      # width.loh <- tapply(q.loh$width, q.loh$loh, FUN=sum)
      # mainq.loh <- names(width.loh)[which.max(width.loh)]
    }
    
    df <- data.frame(cn.p= mainp.cn, 
                     # loh.p = mainp.loh,
                     cn.q= mainq.cn 
                     # loh.q = mainq.loh
                     )
  })
  
  chr.df <- do.call(rbind,chr.list)
  colnames(chr.df) <- paste0(chr,".",colnames(chr.df))
  chr.df <- chr.df[rownames(chrs.df),]
  chrs.df <- cbind(chrs.df,chr.df)
}

# chrs.df <- chrs.df[,grepl("cn",colnames(chrs.df))]
chrs.df <- chrs.df[scoring.table$Sample.names,]
chrs.df$Groups <- scoring.table$Groups_order
# chrs.df <- chrs.df[,grepl("cn",colnames(chrs.df))]
g3 <- chrs.df[chrs.df$Groups ==3,]
barplot(g3)

g3 <- apply(g3,2,function(x){
  ifelse(is.na(x), TRUE,x)})
g3 <- apply(g3,2,function(x){
  ifelse(x > 2, TRUE,FALSE)})
g3 <- g3[,-ncol(g3)]
g3.freq <- apply(g3,2,table)
barplot(g3.freq[c(2,1),])

dev.off()

chrs.df <- apply(chrs.df,2,function(x){
  ifelse(is.na(x), TRUE,x)})
chrs.df <- apply(chrs.df,2,function(x){
  ifelse(x > 2, TRUE,ifelse(x<=2, FALSE,x))})

barplot(g3,beside = T)



scoring.table$Groups_order
chrs.df <- chrs.df[scoring.table$Sample.names,]
chrs.df <- as.data.frame(chrs.df)
chrs.df$Groups <- scoring.table$Groups_order
chrs.df[sapply(chrs.df, is.logical)] <- lapply(chrs.df[sapply(chrs.df,is.logical)],as.factor)

freq.tab <- chrs.df
freq.tab <- freq.tab[,c(ncol(freq.tab),1,2,3,4,5,6,7,8)]
groups <- c("1","2","3","4")
# names(groups) <-c("#DA6220","tomato","goldenrod","#E09B9B")
names(groups) <- c("#7C2C00",  "tomato","goldenrod", "#E09B9B")

freq.group <- data.frame(matrix(nrow = 2,ncol = ncol(freq.tab)-1))
rownames(freq.group) <- c("FALSE", "TRUE")
colnames(freq.group) <- colnames(freq.tab)[-1]
# png(file.path("./images","Genomic_characterization_groups_AllElements.png"),width = 1300,height = 300)
svg(file.path("./images","Genomic_characterization_groups_ChrsSelected.svg"),width = 20,height =4)
par(mfrow=c(1,4), mar=c(13,5,5,5))
j=1
freq.gr.list <- list()
for( j in seq_len(length(groups))){
  grp <- groups[j]
  group.t <- freq.tab[freq.tab$Groups == grp,]
  n.samples <- nrow(group.t)
  group.t <- group.t[,-1]
  i =1
  for(i in seq_len(ncol(group.t))){
    element <- colnames(group.t)[i]
    element.table <- table(group.t[,element])
    freq.group[,element] <- (element.table/n.samples)*100
    freq.group <- freq.group
  }
  freq.gr <- data.matrix(freq.group[c(2,1),])
  freq.gr.list[[paste0("Group_",grp)]]<- freq.gr
  barplot(freq.gr, main = paste0("Group ", grp, " N=",n.samples), col = c(names(grp),"grey90"), border = NA,
          ylim = c(0,100), cex.axis = 2,las =2, cex.names = 1.5,cex.main =2, cex.lab = 2,
          ylab = "%")
}

dev.off()

freq.table <- data.frame(do.call(rbind,freq.gr.list))
freq.table <- apply(freq.table,2,function(x){
  ifelse(x>50,TRUE,FALSE)
})
freq.table <- as.data.frame(freq.table)
freq.table$Groups <-  rep(groups, each =2)
ft <- freq.table[grepl("TRUE",rownames(freq.table)),]
ft <- apply(ft,2,function(x){
  ifelse(x>50,TRUE,FALSE)
})
barplot(data.matrix(freq.table)[,1:8]~Groups,data = data.matrix(freq.table),las=2)

colnames(chrs.df)
df <- chrs.df %>%
  group_by(Groups,chr2.cn.p)%>%
  summarize(Freq =n()) %>%
  data.frame()


df <- chrs.df %>% 
  group_by(Groups) %>%
  summarise(across(starts_with("chr"),
                    list(
                      TRUE_count = ~sum(. == TRUE),
                      FALSE_count = ~sum(. ==FALSE)
                    ),.names = "{.col}_{.fn}"
                    )
            )%>%
  data.frame()

df <- t(data.matrix(df))
colnames(df) <- df[1,]
df <-df[-1,]

        barplot(df,beside  =2)
        

s.g1 <- chr.df

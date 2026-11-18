library(ComplexHeatmap)
library(circlize)
library(dichromat)
library(RColorBrewer)
library(gplots)
library(karyoploteR)
library(dplyr)

score.tot.table <- read.table(file.path("SUPPLEMENTARY_TABLE_S8"), sep="\t", header = TRUE )[,-1]
score.tot.table <- score.tot.table %>%
  mutate(across(c("Purity.est","CNV.freq.2","CNV.freq.3","SV.freq","Number.SNVs","LOH.prop","CN.LOH.REG.prop"), ~ as.numeric(gsub(",", ".", .x))))

# score.tot.table$NF1.PRC2<-FALSE
# score.tot.table$NF1.PRC2[score.tot.table$NF1 ==TRUE & score.tot.table$PRC2 == TRUE] <- TRUE
# score.tot.table$NF1.CDKN2A <- FALSE
# score.tot.table$NF1.CDKN2A[score.tot.table$NF1 ==TRUE & score.tot.table$CDKN2A == TRUE] <- TRUE
# score.tot.table$CDKN2A.PRC2 <- FALSE
# score.tot.table$CDKN2A.PRC2[score.tot.table$PRC2 ==TRUE & score.tot.table$CDKN2A == TRUE] <- TRUE

# View(score.tot.table[,c("TSG.groups","Activating","TSG","NF1.PRC2","NF1.CDKN2A","CDKN2A.PRC2","Groups_order")])
# score.tot.table <- scoring.table
score.tot.table$TSG.groups <- 10
score.tot.table$TSG.groups[score.tot.table$NF1.PRC2==TRUE] <-2
score.tot.table$TSG.groups[score.tot.table$NF1.CDKN2A==TRUE] <-3
score.tot.table$TSG.groups[score.tot.table$CDKN2A.PRC2==TRUE] <-4
score.tot.table$TSG.groups[score.tot.table$TSG==TRUE] <-1

# tt <- score.tot.table[order( !score.tot.table$Activating, score.tot.table$TSG,
#                             score.tot.table$NF1.CDKN2A,score.tot.table$NF1.PRC2, score.tot.table$CDKN2A.PRC2,
#                             score.tot.table$SUZ12, score.tot.table$EED, score.tot.table$TP53, 
#                             score.tot.table$LOH.prop,score.tot.table$CN.LOH.REG.prop),]
tt <- score.tot.table[order(score.tot.table$Activating, score.tot.table$TSG.groups,
                             !score.tot.table$SUZ12, score.tot.table$EED, score.tot.table$TP53, 
                             score.tot.table$LOH.prop,score.tot.table$CN.LOH.REG.prop),]
rownames(tt)<- tt$Sample.names
#NoCellLines
tt <- tt[grepl("MPNST",rownames(tt)),]
#onlyCell lines
# tt <- tt[!grepl("MPNST",rownames(tt)),]

rownames(tt)<- tt$Paper.names

tt$Groups_order[tt$Groups_order ==5] <- 4
tt$Sample.set <- tt$Paper.names
tt$Sample.set[!grepl("Tumor",tt$Sample.set)] <- "GeM"
tt$Sample.set[grepl("Tumor",tt$Sample.set)] <- "Institutional"

# scoring.table$Groups_order
data.test <- tt
data.test <- data.frame(t(data.test))
# data.test <- data.test[c(2:(nrow(data.test))),]

rownames(data.test)
#Prepare the data to plot
#MPNST score
mpnst.score <- as.numeric(data.test["MPNST.score",])
names(mpnst.score) <- colnames(data.test)
mpnst.score.df <- t(data.frame(mpnst.score))
#We delete this column from data test
# data.test <- data.test[-2,]

class(data.test$MPNST.score)
#SV.freq
SV.freq <- data.test["SV.freq",]
SV.freq[is.na(SV.freq)] <- 0
as.numeric(SV.freq)

#CNV.2n.freq
cnv.2n.freq <- data.test["CNV.freq.2",]
cnv.2n.freq[is.na(cnv.2n.freq)] <- 0

#CNV.2n.freq
cnv.3n.freq <- data.test["CNV.freq.3",]
cnv.3n.freq[is.na(cnv.3n.freq)] <- 0

# Number SNVs
nsnvs <- data.test["Number.SNVs",]

#MPNST score annotation barplot
top_annot = HeatmapAnnotation(
  Scoring = anno_barplot(as.numeric(mpnst.score), 
                         bar_width = 1, 
                         gp = gpar(col = "white", fill =c(rep("#D42074",27),rep("grey85",length(28:ncol(mpnst.score.df))))),
                         border = FALSE,
                         axis_param = list(at = c(-10, -5, 0, 5, 10),
                                           labels = c("-10", "-5", "0", "5", "10")),
                         height = unit(2, "cm"), 
                         # col = list(Class = c(rep("#D42074",27),rep("grey85",length(28:ncol(mpnst.score.df)))))
  ),
  # "#.SNVs" = anno_barplot(as.numeric(nsnvs), 
  #                        bar_width = 1, 
  #                        gp = gpar(col = "white", fill =c("#366F6F")),
  #                        border = FALSE,
  #                        axis_param = list(at = c(0, 50000, 100000, 150000, 200000, 250000),
  #                                          labels = c("0", "50000", "100000", "150000", "200000", "250000")),
  #                        height = unit(1, "cm"), 
  #                        # col = list(Class = c(rep("#D42074",27),rep("grey85",length(28:ncol(mpnst.score.df)))))
  # ),
  # SV.frequency = anno_barplot(as.numeric(SV.freq), 
  #                     bar_width = 1, 
  #                     gp = gpar(col = "white", fill =c("#366F6F")),
  #                     border = FALSE,
  #                     axis_param = list(at = c(0,0.5,1),
  #                                       labels = c( "0", "0.5", "1")),
  #                     height = unit(1, "cm"), 
  #                     # col = list(Class = c(rep("#D42074",27),rep("grey85",length(28:ncol(mpnst.score.df)))))
  # ),
  # CNV.2n.freq = anno_barplot(as.numeric(cnv.2n.freq), 
  #                     bar_width = 1, 
  #                     gp = gpar(col = "white", fill =c("#E18B37")),
  #                     border = FALSE,
  #                     axis_param = list(at = c(0,0.5,1),
  #                                       labels = c( "0", "0.5", "1")),
  #                     height = unit(1, "cm"), 
  #                     # col = list(Class = c(rep("#D42074",27),rep("grey85",length(28:ncol(mpnst.score.df)))))
  # ),
  # CNV.3n.freq = anno_barplot(as.numeric(cnv.3n.freq), 
  #                     bar_width = 1, 
  #                     gp = gpar(col = "white", fill =c("#E18B37")),
  #                     border = FALSE,
  #                     axis_param = list(at = c(0,0.5,1),
  #                                       labels = c( "0", "0.5", "1")),
  #                     height = unit(1, "cm"), 
  #                     # col = list(Class = c(rep("#D42074",27),rep("grey85",length(28:ncol(mpnst.score.df)))))
  # ),
  show_annotation_name = T
  # gap =  unit.c(unit(2, "cm"), unit(1, "cm"), unit(1, "cm"))
)
# Specify custom heights for each annotation and add a space

plot(top_annot)


####### Bottom annotation
#purity degradation annotation
# blues9
# rownames(data.test)
# col_fun = colorRamp2(c(0, 25, 50, 75, 100), c(viridis::viridis(5)))
col_fun = colorRamp2(c(50, 75, 100), c("#F7FBFF","#6DAFD6","#08306B"))
Purity.est <- as.numeric(data.test["Purity.est",])
col_fun_2 = colorRamp2(c(0, 0.3, 0.65), c("#FAFCFC","#8CACAC","#366F6F"))
# max(as.numeric(CNV.freq.2),na.rm = T)

col_fun_3 = colorRamp2(c(0, 0.5, 1),c("#FDFBF6","#EDC88C","#E18B37"))
CNV.freq.2 = as.numeric(data.test["CNV.freq.2",])
col_fun_4 = colorRamp2(c(0, 0.5, 1), c("#FDFBF6","#EDC88C","#E18B37"))
CNV.freq.3 = as.numeric(data.test["CNV.freq.3",])

col_fun.5 = colorRamp2(c(0,100000,200000),c("#E6EBFE","#3A56B3","#0B1123") )
SNV.num = as.numeric(data.test["Number.SNVs",])

# RColorBrewer::brewer.pal(5,"RdYlBu")
# RColorBrewer::brewer.pal.info
#LOH.prop
cn.loh.reg.prop <- as.numeric(data.test["CN.LOH.REG.prop",])
cn.loh.reg.prop <- data.frame(cn.loh.reg.prop)
rownames(cn.loh.reg.prop) <- colnames(data.test)
cn.loh.reg.prop$cn.loh.reg.prop[cn.loh.reg.prop$cn.loh.reg.prop<=50] <- 0
cn.loh.reg.prop$cn.loh.reg.prop[cn.loh.reg.prop$cn.loh.reg.prop>50] <- 1

loh.prop <- as.numeric(data.test["LOH.prop",])
loh.prop <- data.frame(loh.prop)
rownames(loh.prop) <- colnames(data.test)
loh.prop$loh.prop[loh.prop$loh.prop<=50] <- 0
loh.prop$loh.prop[loh.prop$loh.prop>50] <- 1

bottom_annot <- cbind(cn.loh.reg.prop,loh.prop)
rownames(data.test)
# data.test.2 <- data.test[c(15,1,16:18,19:23,26:28,31:39),]
# data.test.2 <- data.test[c(15,26,20:24,27:29,32),]
# data.test.2 <- data.test[c(19:23,26:28,31:39,1),]
# data.test.2 <- data.test[c("NF1","CDKN2A","CDKN2A.SV","SUZ12","EED","TP53","RB1","Activating",
#                            "chr1.unbalance","cn.loh.reg.prop","loh.prop","Main.Dx.Histo", "Histology.pattern","Methylation","Purity.est","Type"),]
data.test.2 <- data.test[c("Groups_order","Cluster", "Sample.set", "Sex","TSG","NF1","CDKN2A","CDKN2A.SV","SUZ12","EED","TP53","RB1","PTEN","Activating",
                           "chr1.unbalance","Main.Dx.Histo", "Histology.pattern","Methylation","H3K27me3","Purity.est","Type","Number.SNVs","SV.freq","CNV.freq.2","CNV.freq.3"),]

# groups <- c("#DA6220","tomato","goldenrod","#C55656" ,"#E09B9B")
# names(groups) <- c("1","2","3","4","5")
# groups <- c("#DA6220","tomato","goldenrod","#E09B9B")
groups <- c("#7C2C00",  "tomato","goldenrod", "#E09B9B")
names(groups) <- c("1","2","3","4")
colnames(data.test.2)
# bottom_annot[18]
data.test.2 <- t(data.test.2)
bottom_annot <- cbind(data.test.2,cn.loh.reg.prop,loh.prop)
# bottom_annot <- cbind(bottom_annot,t(data.test[c(16:19),]), Purity.est=as.numeric(Purity.est), t(data.test[c(1),]))
bottom_annot <- bottom_annot[,c("Groups_order", "Cluster", "Sample.set", "Type", "Sex", "TSG","NF1","CDKN2A","CDKN2A.SV","SUZ12","EED","TP53","RB1","PTEN","Activating",
                                "chr1.unbalance","cn.loh.reg.prop","loh.prop", "SV.freq","CNV.freq.2","CNV.freq.3","Number.SNVs","Main.Dx.Histo", "Histology.pattern","H3K27me3","Methylation","Purity.est")]
bottom_annot$Purity.est <- as.numeric(bottom_annot$Purity.est)
bottom_annot$SV.freq<- as.numeric(bottom_annot$SV.freq)
bottom_annot$CNV.freq.2<- as.numeric(bottom_annot$CNV.freq.2)
bottom_annot$CNV.freq.3<- as.numeric(bottom_annot$CNV.freq.3)
bottom_annot$Number.SNVs<- as.numeric(bottom_annot$Number.SNVs)


color.truefalse <- c( "TRUE"="#D8BDF9","FALSE" = "grey95")

mycolors <- list(Type = c("NF1" ="#5AC1D1","SP"="#815AD1"),
                 Sex = c("F" = "#B3AD3A","M" ="#3A40B3"),
                 Groups_order =groups,
                 Sample.set = c("Institutional"="#0B007A", GeM ="#87CEFA" ),
                 Cluster = c(C1 = "#5279FB","Grey"="grey55" , C2 = "#FBD352"),
                 # Class=c(MPNST = "#D42074", MissClass="grey85"),
                 NF1=color.truefalse,
                 CDKN2A=color.truefalse,
                 CDKN2A.SV = c("TRUE"="#8E3DEE","FALSE"="grey95"),
                 SUZ12=color.truefalse,
                 EED=color.truefalse,
                 PRC2=color.truefalse,
                 TSG= c( "TRUE"="#FF6625","FALSE" = "grey95"),
                 TP53=color.truefalse,
                 RB1=color.truefalse,
                 PTEN =color.truefalse,
                 Activating=c("TRUE"="#0B007A","FALSE"="grey95"),
                 chr1.unbalance = c("TRUE"= "#87CEFA","FALSE"="grey95"),
                 # chr1.p.cn = c("cn.loh"= "#87CEFA", "gain"= "#F49738","het.loss" = "#65FF8F","2n" = "#E5E5E5" ),
                 # chr1.q.cn = c("cn.loh"= "#87CEFA", "gain"= "#F49738","het.loss" = "#65FF8F","2n" = "#E5E5E5" ),
                 # chr8.unbalance = color.truefalse,
                 # chr8.p.cn = c("cn.loh"= "#87CEFA", "gain"= "#F49738","het.loss" = "#65FF8F","2n" = "#E5E5E5" ),
                 # chr8.q.cn = c("cn.loh"= "#87CEFA", "gain"= "#F49738","het.loss" = "#65FF8F","2n" = "#E5E5E5" ),
                 # chr17.unbalance = color.truefalse,
                 # chr17.p.cn = c("cn.loh"= "#87CEFA", "gain"= "#F49738","het.loss" = "#65FF8F","2n" = "#E5E5E5" ),
                 # chr17.q.cn = c("cn.loh"= "#87CEFA", "gain"= "#F49738","het.loss" = "#65FF8F","2n" = "#E5E5E5" ),
                 # 
                 Main.Dx.Histo = c("04=MPNST-low grade"="#52FBB1",
                                   "05=Conventional high grade MPNST (hyper- and hypo-cellular fascicles with/without geographic necrosis)" = "#FB529C",
                                   "06=Conventional MPNST with heterologous elements" = "#FBB152",
                                   "07=Non-conventional MPNST with associated conventional low grade areas" = "#B052FB",
                                   "08=Non-conventional MPNST"= "#529CFB",
                                   "09=Spindle cell/undifferentiated sarcoma"="#101423"),
                 Histology.pattern = c("01=Rhabdomyoblastic (Triton)"= "#529CFB",
                                       "01=Rhabdomyoblastic (Triton); 02=Osteo/chodrosarcomatous"="#87CEFA",
                                       "02=Osteo/chodrosarcomatous"="#FBB152",
                                       "03=Angiosarcomatous" = "#DC143C",
                                       "Not applicable" ="grey"),
                 Methylation = c(neurofibroma = "#52FBB1",
                                 "epitheliod sarcoma"="#FFE4E1",
                                 "MPNST" = "#B4077E",
                                 "MPNST-Like" = "#DE94C7",
                                 "osteosarcoma" ="#FBB152",
                                 "rhabdomyosarcoma" ="#529CFB",
                                 "Undifferentiated Sarcoma" = "black",
                                 "Well differentiated liposarcoma" = "#FFE821",
                                 "Melanoma" = "#784A4D",
                                 "NS" = "grey"),
                 H3K27me3 = c(Loss="#52B5FB", Retained="#FB5272"),
                 Purity.est = col_fun,
                 SV.freq = col_fun_2,
                 CNV.freq.2 = col_fun_3,
                 CNV.freq.3 = col_fun_4,
                 Number.SNVs = col_fun.5,
                 cn.loh.reg.prop=c("1"="#CD853F","0"="grey95"), 
                 loh.prop=c("1"="#6591A5","0"="grey95"))

ha2 = HeatmapAnnotation(df =bottom_annot,
                        col = mycolors,gp = gpar(col = "white"))


plot(ha2)
dev.off()


bottom_annot$Purity.est



#Scoring values

scoring.val <- data.test[grepl("^s.",rownames(data.test)),]
scoring.val.df <- apply(scoring.val,2,as.numeric)
rownames(scoring.val.df) <- rownames(scoring.val)
scoring.val.df <- scoring.val.df[-c((nrow(scoring.val.df)-1),nrow(scoring.val.df)),]
scoring.val.df[1,] <- 0
scoring.val.df <- t(data.frame(scoring.val.df[1,]))

#scoring colors
scoring.col = colorRamp2(c(-3,0,3), c( "cornflowerblue","white", "#FA9F57"))
# annot.scoring <- data.frame(scorign.values = rownames(scoring.val.df), groups = c(rep("TSG",5),"Activating","LOH","LOH","CNV"))
# rownames(annot.scoring) <-annot.scoring$scorign.values
# annot.scoring <- annot.scoring$groups
# names(annot.scoring) <- rownames(scoring.val.df)
# annot.scoring <- data.frame(groups = factor(annot.scoring, levels=c("TSG","Activating","LOH","CNV")))
# annot.col <- list(groups = c(TSG = "#BF5CFF",Activating = "#FF6625", LOH = "cornflowerblue", CNV = "#FFE394"))
# 
# ha3 = rowAnnotation(df =annot.scoring,
#                     col = annot.col)
# plot(ha3)

heat.map.plot <- Heatmap(matrix = scoring.val.df, row_title = "",name = "Scoring",
                         cluster_columns = F,
                         cluster_rows = F,
                         width = ncol(scoring.val.df)*unit(5, "mm"), 
                         height = nrow(scoring.val.df)*unit(5, "mm"),
                         rect_gp = gpar(col= "white"),
                         col = scoring.col,
                         # top_annotation = top_annot,
                         bottom_annotation = ha2,
                         # left_annotation = ha3,
                         column_gap = unit(3, "mm"))

svg(filename = file.path("./images/Paper.Names.table.v13_Order_test_v15.svg"), width = 25, height = 20)
draw(heat.map.plot, ht_gap = unit(3, "mm"))
dev.off()

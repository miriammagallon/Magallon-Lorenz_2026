###############################################
# Chi square Genomic Groups characterization  #
###############################################
library(dplyr)
#Packages needed

#loading information

 T)
scoring.table <- read.table("Supplementary_Table_S8", sep ="\t",header = T)
scoring.table <- scoring.table %>%
  mutate(across(c("Purity.est","CNV.freq.2","CNV.freq.3","SV.freq","Number.SNVs","LOH.prop","CN.LOH.REG.prop"), ~ as.numeric(gsub(",", ".", .x))))

#NoMPNST cell lines
scoring.table <- scoring.table[grepl("MPNST",scoring.table$Sample.names),]


# scoring.table$Groups <-"0" 
# scoring.table$Groups[1:38] <- "1"
# scoring.table$Groups[39:46] <- "0"
# scoring.table$Groups[47:53] <- "2"
# scoring.table$Groups[54:59] <- "3"
# colnames(scoring.table)
freq.tab <- data.frame(Sample.name=scoring.table$Sample.name, Group= scoring.table$Groups_order, EED= scoring.table$EED, Activating =scoring.table$Activating,
                       "TSG 3KO" = scoring.table$TSG, TP53 =scoring.table$TP53, CDKN2A = scoring.table$CDKN2A, chr8.unbal = scoring.table$chr8.unbalance,
                       PRC2 = scoring.table$PRC2, NF1 = scoring.table$NF1, CDKN2A.SV = scoring.table$CDKN2A.SV, SUZ12 =scoring.table$SUZ12, Cluster = scoring.table$Cluster,
                       "LOH" = scoring.table$LOH.prop, "CN-LOH regions" = scoring.table$CN.LOH.REG.prop, "chr1 unbal"=scoring.table$chr1.unbalance, Sex = scoring.table$Sex)
rownames(freq.tab) <- freq.tab$Sample.name
freq.tab <- freq.tab[,-1]

# freq.tab$group <- factor(freq.tab$group, levels = c("1","0","2","3","4"))
freq.tab[sapply(freq.tab, is.logical)] <- lapply(freq.tab[sapply(freq.tab,is.logical)],as.factor)
# freq.tab$cd

# #Inactivation of activating in G1 vs G0
colnames(freq.tab)
# freq.tab <- freq.tab[,c(1,12,4,9,6,10,8,11,2,5,3,15,13,14)]
freq.tab <- freq.tab[,c(1,12,4,9,6,10,5,8,11,2,3)] #without the genomic structure part



#build a frequency table with all elements of group.t
# groups <- c("1","2","3","0")
# names(groups) <- c("#DA6220",  "tomato","goldenrod", "#E09B9B")

groups <- c("1","2","3","4")
# names(groups) <-c("#DA6220","tomato","goldenrod","#E09B9B")
names(groups) <- c("#7C2C00",  "tomato","goldenrod", "#E09B9B")

freq.tab$Group[freq.tab$Group == 5] <- 4
freq.group <- data.frame(matrix(nrow = 2,ncol = ncol(freq.tab)-2))
rownames(freq.group) <- c("FALSE", "TRUE")
colnames(freq.group) <- colnames(freq.tab)[-c(1,2)]

svg(file.path("./images","Genomic_characterization_groups_AllElements_v8.svg"),width = 20,height =4)
par(mfrow=c(1,4), mar=c(13,5,5,5))
j=5
for( j in seq_len(length(groups))){
  grp <- groups[j]
  group.t <- freq.tab[freq.tab$Group == grp,-c(2)]
  n.samples <- nrow(group.t)
  i =12
  for(i in seq_len(ncol(group.t)-1)){
    element <- colnames(group.t)[i+1]
    element.table <- table(group.t[,element])
    freq.group[,element] <- (element.table/n.samples)*100
    freq.group <- freq.group
  }
  freq.gr <- data.matrix(freq.group[c(2,1),])
  barplot(freq.gr, main = paste0("Group ", grp, " N=",n.samples), col = c(names(grp),"grey90"), border = NA,
          ylim = c(0,100), cex.axis = 2,las =2, cex.names = 1.5,cex.main =2, cex.lab = 2,
          ylab = "%")
}
 dev.off()






##################################
# Pie charts histology by group  #
##################################

library(dplyr)
groups <- c("1","2","3","4")
names(groups) <- c("#7C2C00",  "tomato","goldenrod", "#E09B9B")
# scoring.table$Groups_order[scoring.table$Groups_order==5] <- 4 
levels.df <- unique(scoring.table$Main.Dx.Histo)
levels.df <-levels.df[c(2,1,4:7,3)]
df.freq <- data.frame(Groups =scoring.table$Groups_order, 
                      Sex = scoring.table$Sex,
                      Main.Hist = scoring.table$Main.Dx.Histo,
                      NF1vsSP =scoring.table$Type,
                      chr1.unbalance = scoring.table$chr1.unbalance, 
                      Cluster = scoring.table$Cluster,
                      PRC2 = scoring.table$PRC2,
                      SUZ12 = scoring.table$SUZ12,
                      EED = scoring.table$EED,
                      TP53 = scoring.table$TP53,
                      CDKN2A = scoring.table$CDKN2A,
                      CDKN2A.SV = scoring.table$CDKN2A.SV,
                      SV.freq = scoring.table$SV.freq,
                      LOH.prop = scoring.table$LOH.prop,
                      CNV.freq.2 = scoring.table$CNV.freq.2,
                      Gain.freq = scoring.table$CNV.freq.3,
                      Methilome.Classifier = scoring.table$Methylation)
df.freq$Cluster <- factor(df.freq$Cluster, levels = c("C1", "C2", "Grey"))
# df.freq$Main.Hist[is.na(df.freq$Main.Hist)] <- "NA"
df.freq$Main.Hist <- factor(df.freq$Main.Hist, levels = c(levels.df))
df.freq$Sex <- factor(df.freq$Sex, levels = c("M","F"))
df.freq[sapply(df.freq, is.logical)] <- lapply(df.freq[sapply(df.freq,is.logical)],as.factor)

df <- df.freq %>%
  group_by(Groups,Main.Hist)%>%
  summarize(Freq =n()) %>%
  data.frame()

df$Groups
Main.Dx.Histo = c("05=Conventional high grade MPNST (hyper- and hypo-cellular fascicles with/without geographic necrosis)" = "#FB529C",
                  "04=MPNST-low grade"="#52FBB1",
                  "06=Conventional MPNST with heterologous elements" = "#FBB152",
                  "07=Non-conventional MPNST with associated conventional low grade areas" = "#B052FB",
                  "08=Non-conventional MPNST"= "#529CFB",
                  "09=Spindle cell/undifferentiated sarcoma"="#101423"
                  )
NFvsSP <- c("NF1" ="#5AC1D1","SP"="#815AD1")

svg(file.path("./images","PieCharts_histologyBygroups_v7.svg"),width = 40,height = 15)
par(mfrow=c(1,5), mar=c(5,5,5,5))
i=1
for(i in seq_len(length(groups))){
  g <- groups[i]
  df.1 <- df[df$Groups == g,]
  data.summary <- aggregate(Freq~ Main.Hist, df.1,sum)
  
  colours.df <- Main.Dx.Histo[names(Main.Dx.Histo)%in%data.summary$Main.Hist]
  
  pie(data.summary$Freq, col = colours.df,
      main = paste0("Group ",g, " N=",sum(data.summary$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
  
}
par(mfrow=c(1,1), mar=c(5,5,5,5))

legend("bottom",ncol = 1,bty = "n",border = NA, legend = names(Main.Dx.Histo), fill = Main.Dx.Histo, cex =3)

dev.off()

# 
# 
# df.freq <- data.frame(Groups =scoring.table$Groups_order, Cluster = scoring.table$Cluster)
# df.freq$Cluster <- factor(df.freq$Cluster, levels = c("C1", "C2", "Grey")), Cluster = scoring.table$Cluster)
# df.freq$Cluster <- factor(df.freq$Cluster, levels = c("C1", "C2", "Grey"))

df <- df.freq %>%
  group_by(Groups,Cluster)%>%
  summarize(Freq =n()) %>%
  data.frame()


Cluster = c(C1 = "#5279FB", C2 = "#FBD352","Grey"="grey55" )

svg(file.path("./images","PieCharts_ClusterByGroups_v7.svg"),width = 40,height = 15)
par(mfrow=c(1,5), mar=c(5,5,5,5))
i=1
for(i in seq_len(length(groups))){
  g <- groups[i]
  df.1 <- df[df$Groups == g,]
  data.summary <- aggregate(Freq~ Cluster, df.1,sum)
  
  colours.df <- Cluster[data.summary$Cluster]
  
  pie(data.summary$Freq,col = colours.df,
      main = paste0("Group ",g, " N=",sum(data.summary$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
  
}
par(mfrow=c(1,1), mar=c(5,5,5,5))

legend("bottom",ncol = 1,bty = "n",border = NA, legend = names(Cluster), fill = Cluster, cex =3)

dev.off()

#Sex

df <- df.freq %>%
  group_by(Groups, Sex)%>%
  summarize(Freq =n()) %>%
  data.frame()


Sex = c("F" = "#B3AD3A","M" ="#3A40B3", "NA" = "grey")

svg(file.path("./images","PieCharts_SexByGroups_v7.svg"),width = 40,height = 15)
par(mfrow=c(1,5), mar=c(5,5,5,5))
i=1
for(i in seq_len(length(groups))){
  g <- groups[i]
  df.1 <- df[df$Groups == g,]
  data.summary <- aggregate(Freq~ Sex, df.1, sum)
  
  colours.df <- Sex[data.summary$Sex]
  
  pie(data.summary$Freq,col = colours.df,
      main = paste0("Group ",g, " N=",sum(data.summary$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
  
}
par(mfrow=c(1,1), mar=c(5,5,5,5))

legend("bottom",ncol = 1,bty = "n",border = NA, legend = names(Sex), fill = Sex, cex =3)

dev.off()



df <- df.freq %>%
  group_by(Groups,NF1vsSP)%>%
  summarize(Freq =n()) %>%
  data.frame()
svg(file.path("./images","PieCharts_NFvsSPByGroups_v7.svg"),width = 40,height = 15)
par(mfrow=c(1,5), mar=c(5,5,5,5))
i=3
for(i in seq_len(length(groups))){
  g <- groups[i]
  df.1 <- df[df$Groups == g,]
  if(nrow(df.1)==1){
    pie(df.1$Freq,col = colours.df,
        main = paste0("Group ",g, " N=",sum(df.1$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
    
  }else{
    data.summary <- aggregate(Freq~ NFvsSP, df.1,sum)
    
    colours.df <- NFvsSP
    
    pie(data.summary$Freq,col = colours.df,
        main = paste0("Group ",g, " N=",sum(data.summary$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
    
  }
 
}
par(mfrow=c(1,1), mar=c(5,5,5,5))

legend("bottom",ncol = 1,bty = "n",border = NA, legend = names(NFvsSP), fill = NFvsSP, cex =3)

dev.off()                     


df <- df.freq %>%
  group_by(Groups,chr1.unbalance)%>%
  summarize(Freq =n()) %>%
  data.frame()
df$chr1.unbalance <- factor(df$chr1.unbalance,levels=c("TRUE","FALSE"))
unbalance = c("TRUE" = "#97667B", "FALSE" = "grey85" )


svg(file.path("./images","PieCharts_chr1.UnbalanceByGroups_v7.svg"),width = 40,height = 15)
par(mfrow=c(1,5), mar=c(5,5,5,5))
i=3
for(i in seq_len(length(groups))){
  g <- groups[i]
  df.1 <- df[df$Groups == g,]
  if(nrow(df.1)==1){
    pie(df.1$Freq,col = colours.df,
        main = paste0("Group ",g, " N=",sum(df.1$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
    
  }else{
    data.summary <- aggregate(Freq~ chr1.unbalance, df.1,sum)
    
    colours.df <- unbalance
    
    pie(data.summary$Freq,col = colours.df,
        main = paste0("Group ",g, " N=",sum(data.summary$Freq)),labels = NA, cex.main = 6, clockwise = TRUE, init.angle = 90)
    
  }
  
}
par(mfrow=c(1,1), mar=c(5,5,5,5))

legend("bottom",ncol = 1,bty = "n",border = NA, legend = names(unbalance), fill = unbalance, cex =3)

dev.off()    
df.freq
# svg(filename=file.path(execution.dir,"results","images", paste0("BoxPlot_CNVprop_cluster1.2_v3.svg")), width = 7, height = 8)

tabla <- table(df.freq$chr1.unbalance, df.freq$Groups)
tabla <- tabla[c(2,1),]
svg(file.path("./images","barplot_chr1.UnbalanceByGroups_v7.svg"),width = 7, height = 8)
# par(mfrow=c(1,1), mar=c(5,5,5,5))
barplot(tabla,beside = FALSE,col = unbalance, main = "chr1 unbalance", cex.axis = 3, cex.names = 3,horiz = F, ylim = c(0,35),las=2, names.arg = c("G1","G2","G3","Others"))
legend("topright",ncol = 1,bty = "n",border = NA, legend = names(unbalance), fill = unbalance, cex =3)
dev.off()    


##### Statistical test review. Fisher test

# Función para hacer test de Fisher por pares en cualquier tabla de contingencia 2D
pairwise_fisher_base <- function(mat, p.adjust.method = "bonferroni") {
  col_names <- colnames(mat)
  # message( "pairwise fisher ", col_names)
  k <- ncol(mat)
  p_vals <- c()
  pairs <- c()
  
  for(i in 1:(k-1)) {
    for(j in (i+1):k) {
      # Extrae la subtabla para cada par de grupos (columnas)
      sub_mat <- mat[, c(i, j)]
      # Ejecuta el test de Fisher (simulated.p.value es útil si hay casillas con ceros)
      test <- fisher.test(sub_mat, simulate.p.value = TRUE) 
      p_vals <- c(p_vals, test$p.value)
      pairs <- c(pairs, paste("Group", col_names[i], "vs", col_names[j]))
    }
  }
  
  # Aplica la corrección de p-valores (ej. Bonferroni)
  p_adj <- p.adjust(p_vals, method = p.adjust.method)
  
  # Devuelve un data.frame limpio con los resultados
  data.frame(Comparison = pairs, p.raw = p_vals, p.adj = p_adj)
}


totest <- colnames(df.freq)[-c(1,12,13,14,15,16)]
i=2
for(i in seq_len(length(totest))){
  cn <- totest[i]

  tc <- table(df.freq[,cn], df.freq$Groups)
  mat =tc
  print(tc)
  message(cn ," fisher test")
  pfb <- pairwise_fisher_base(tc, p.adjust.method = "bonferroni")
  print(pfb)
}
fisher.test.tc <- fisher.test(tc)
print(fisher.test.tc)
library(rstatix)
f
pairwise_fisher_test(tc[c(1,4),], p.adjust.method = "bonferroni")
chisq_test(tc,correct = T)
pairwise.t.test(x = df.freq$SV.freq, g = df.freq$Groups, p.adjust.method = "bonferroni")
pairwise.t.test(x = df.freq$LOH.prop, g = df.freq$Groups, p.adjust.method = "bonferroni")
pairwise.t.test(x = df.freq$CNV.freq.2, g = df.freq$Groups, p.adjust.method = "bonferroni")
pairwise.t.test(x = df.freq$Gain.freq, g = df.freq$Groups, p.adjust.method = "bonferroni")
pairwise
scoring.table$
i=5
for(i in seq_len(length(totest))){
  cn <- totest[i]
  mpsnt.prc <- df.freq[df.freq[,cn] %in% c("MPNST","MPNST-Like"),]
  tc <- table(df.freq[,cn], df.freq$PRC2)
  message(cn ," fisher test")
  pfb <- pairwise_fisher_base(mat = tc, p.adjust.method = "bonferroni")
  print(pfb)
}
totest
cn <- "Methilome.Classifier"
main.classic <- df.freq[,cn]
main.classic.df <- df.freq
to.false <- which(!grepl("MPNST",main.classic))
to.false <- which("MPNST-Like"!=main.classic)

main.classic.df[,cn]<-as.character(main.classic.df[,cn])
main.classic.df[,cn][to.false]<-FALSE
main.classic.df[,cn]<-as.factor(main.classic.df[,cn])

tc <- table(main.classic.df[,cn], main.classic.df$Groups)
mat =tc
print(tc)
message(cn ," fisher test")
pfb <- pairwise_fisher_base(tc, p.adjust.method = "bonferroni")
# 
# 
# 
# 
# barplot(freq.group, main = main, col = c("#97667B","grey90"), border = NA,
#         ylim = c(0,30), cex.axis = 2,las =2, cex.names = 2,cex.main =2, cex.lab = 1.5,
#         ylab = "# Samples")
# 
# par(mfrow=c(2,5))
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,i+1)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "CDKN2A", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,10)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "CDKN2A.SV", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,8)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "PRC2", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,4)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "TSG", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,2)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "EED", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,5)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "TP53", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,3)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "activating mutations", col = c("#97667B","grey90"), border = NA)
# 
# freqtt <- t(table(freq.tab[freq.tab$group%in%c("1","0","2","3","4"),c(1,7)]))
# freqtt <- freqtt[c(2,1),c(2,3,4,5,1)]
# barplot(freqtt,main = "chr8.unbalance", col = c("#97667B","grey90"), border = NA)
# 
# # dev.off()

# 
# 
# ##### Comparisons global characteristics all samples #####
# par(mfrow=c(2,3))
# #Inactivation of CDKN2A vs TP53
# freqtt <- table(freq.tab[,c(6,5)])
# freqtt
# chis <- chisq.test(freqtt)
# chis
# plot(freqtt, col = c("grey90","#97667B"), main = "TP53 inactivation vs CDKN2A inactivation All samples")
# 
# #Inactivation of TSG vs activating
# freqtt <- table(freq.tab[,c(3,4)])
# freqtt
# chis <- chisq.test(freqtt,)
# chis
# plot(freqtt, col = c("grey90","#97667B"), main = "TSG inactivation vs Activating mutations All samples")
# 
# #Inactivation of PRC2 vs activating
# freqtt <- table(freq.tab[,c(3,8)])
# freqtt
# chis <- chisq.test(freqtt)
# chis
# plot(freqtt, col = c("grey90","#97667B"), main = "PRC2 inactivation vs Activating mutations All samples")
# dev.off()
# ##### Comparisons by group of the different genomics characteristics all samples #####
# 
# #Inactivation of EED in G1 vs G3 
# freqtt <- table(freq.tab[freq.tab$group%in%c("1","3"),c(1,2)])
# freqtt
# chis <- chisq.test(freqtt)
# chis
# plot(freqtt, col = c("grey90","#97667B"), main = "EED inactivation in G3 vs G1")
# 
# #Inactivation of TP53 in G1 vs G3 
# freqtt <- table(freq.tab[freq.tab$Group%in%c("1","3"),c(1,5)])
# freqtt
# chis <- chisq.test(freqtt)
# chis
# plot(freqtt, col = c("grey90","#97667B"), main = "TP53 inactivation in G3 vs G1")
# plot.new()
# # dev.off()
# ###
# # #Inactivation of EED in G1 vs G3
# freqtt <- table(freq.tab[freq.tab$group%in%c("1","4"),c(1,4)])
# freqtt
# chis <- chisq.test(freqtt)
# chis
# plot(freqtt, col = c("#FB6E52","#52FBB1"), main = "EED inactivation in G4 vs G1")

##################### Mut Signatures #############
library(mutSignatures)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
devtools::install_github("dami82/mutSignatures", force = TRUE, build_vignettes = TRUE)
# BiocManager::install("mutSignatures")
library(mutSignatures)
################################################################
#Get 30 different colors
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- c(col_vector)
execution.dir<- "."


sample.data <- read.table(file = file.path(execution.dir,"Sample.info.csv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE,comment.char = "")
sample.data <- sample.data[-c(1,30:32),]
sample.names <- sample.data$Paper.Names
names(sample.names) <- sample.data$Sample.name

mut.dir <-  file.path(execution.dir,"results/MutSignatures")
if(!file.exists(mut.dir))dir.create(mut.dir)

image.dir <- file.path(execution.dir,"results/MutSignatures/images")
if(!file.exists(image.dir))dir.create(image.dir)
# mut.signatures.df <- read.table(file = file.path(execution.dir,"results",
                                                 # paste0("AllSamples_MutSignVariants.v4.csv")),sep = "\t", header=T,stringsAsFactors = F)
stk.var<- c("germ.variants","som.variants")
i=2
mut.signatures.df <- read.table(file = file.path(execution.dir,"results",
                                                 paste0("AllSamples_",stk.var[i],"_MutSignVariants.v2.csv")),sep = "\t", header=T,stringsAsFactors = F)


sample.id  <- unique(mut.signatures.df$SAMPLEID)
sid <-sample.names [names(sample.names)%in% sample.id ]
sid <- sid[sid!=""]

mut.signatures.df <- mut.signatures.df[mut.signatures.df$SAMPLEID%in% names(sid),]
sample.id  <- unique(mut.signatures.df$SAMPLEID)
for(i in seq_len(length(sid))){
  sn <- names(sid)[i]
  mut.signatures.df$SAMPLEID[mut.signatures.df$SAMPLEID%in% sn] <- sid[sn]
}


# mut.signatures.df <- mut.signatures.df[mut.signatures.df$SAMLEID != "STS-26T",]
# 
# mut.signatures.df <- read.table(file = file.path(execution.dir,"results", 
# 
# paste0("AllSamples_MutSignVariants.TMB.csv")),sep = "\t", header=T,stringsAsFactors = F)

x <- filterSNV(dataSet = mut.signatures.df, seq_colNames = c("REF","ALT"))
genome <- BSgenome.Hsapiens.UCSC.hg38
head(x)
# x<- output.data

x <- mutSignatures::attachContext(mutData = x,
                   chr_colName = "CHROM",
                   start_colName = "POS.START",
                   end_colName = "POS.END",
                   # context_colName = "REF",
                   nucl_contextN = 3,
                   BSGenomeDb = genome)

# Remove mismatches
x <- removeMismatchMut(mutData = x,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")
# Compute mutType
x <- attachMutType(mutData = x,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context")

# Mutation Counts
ss <-x
x <- x[!grepl("N",x$mutType),]
cells.signature.counts <- countMutTypes(mutTable = x,
                             mutType_colName = "mutType",
                             sample_colName = "SAMPLEID")

# save(cells.signature.counts, file = file.path(execution.dir,"results/MutSignatures", "MPNSTtumors_MutSignatureCounts_som.variants.v2.RData"))
load(file.path(execution.dir,"results/MutSignatures", "MPNSTtumors_MutSignatureCounts_germ.variants.v2.RData"))
# load(file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts.2.RData"))
# save(cells.signature.counts, file = file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts_NO_26T.RData"))

cells.signature.counts
df<- mutSignatures::as.data.frame(cells.signature.counts)
colSums(df)
#Get cosmic signatures 
cosmic <-mutSignatures::getCosmicSignatures("./COSMIC_v3.3_SBS_GRCh38.txt",forceUseMirror = FALSE )
cosmic <-read.table("./Annotations/COSMIC_v3.3_SBS_GRCh38.txt", sep = "\t", header = T,row.names = "Type")
head(cosmic)
# cosmic.selected <- c("COSMIC.16","COSMIC.7","COSMIC.30","COSMIC.12","COSMIC.19","COSMIC.5")
# cosmic.selected <- c("COSMIC.1","COSMIC.3","COSMIC.5","COSMIC.8","COSMIC.9",
#                      "COSMIC.12","COSMIC.16","COSMIC.19","COSMIC.20",
#                      "COSMIC.25","COSMIC.26","COSMIC.30")
# cosmic.selected <- c("COSMIC.16","COSMIC.7","COSMIC.5")
cosmic.selected <- c("COSMIC.1","COSMIC.6","COSMIC.7","COSMIC.11","COSMIC.30")
#
# # cosmic.selected <- rownames(msig1$distanceMatrix)[rowSums(msig1$distanceMatrix)<0.6]
cosmic.df <- mutSignatures::as.data.frame(cosmic)
cosmic.df <- cosmic.df[,cosmic.selected]
cosmic <- mutSignatures::as.mutation.signatures(cosmic)
#

# how many signatures should we extract?
num.sign <- 5
load(file.path(execution.dir,"results","MutSignatures",
               paste0("MPNSTcells_MutSignatureAnalysis_1000Iter_",num.sign,"sign.RData")))

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
params <-
  mutSignatures::setMutClusterParams(
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 500,               # bootstrapping: usually 500-1000
    num_parallelCores = 4)                # total num of cores to use (parallelization)

# Extract new signatures - may take a while
cell.signature.analysis <-
  decipherMutationalProcesses(input = cells.signature.counts,
                              params = params)

# # # save(cells.signature.counts, file = file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts.RData"))
save(cell.signature.analysis, file = file.path(mut.dir, "MPNST_IB_MutSignatureCounts_ExtractSignatures.RData"))

# Retrieve signatures (results)
cell.signature.sig <- cell.signature.analysis$Results$signatures

# Retrieve exposures (results)
cell.signature.exp <- cell.signature.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(cell.signature.sig, signature =1 , ylim = c(0, 0.05))


png(file.path(image.dir,paste0("de-novosiganturesCellLinesTest",num.sign,".png")),width = 1000,height = 500)
msigPlot(cell.signature.exp, main = "cell.lines")+
  scale_fill_manual(values = c(col_vector[1:80]))

dev.off()


#Another way to represent the cell.signature.exp results
cse <- mutSignatures::as.data.frame(cell.signature.exp)
colnames(cse)
for(i in seq_len(length(sample.names))){
  colnames(cse)[colnames(cse) %in% names(sample.names)[i]] <- sample.names[i] 
} 
#tumor wsg
if(stk.var[i]== "germ.variants"){
  # cse <- cse[,c(8,6,10,13,15,17,19,21,25,28,1,23,2,3,4,5,7,9,11,12,14,16,18,20,22,24,26,27,29)] 
  cse <- cse[,c(7,9,12,14,16,26,1,2,6,3,4,5,34,11,13,15,17,18,19,21,23,24,25,27,30,32,33,35,20,8,10,22,29,28)]
  # cse <- cse[,c(1,2,3,4,34,6,5,17,18,19:21,23,24,25,27,30,32,33,35,36,31,22,29,28)] #IB selected samples
}

df.sig <- c()
for(i in seq_len(ncol(cse))){
  sn <- colnames(cse)[i]
  df.sig<- c(df.sig,cse[,sn])
  
}
df.sig <- data.frame("Mutation counts by signature" = df.sig, Signatures = rep(rownames(cse),ncol(cse)))
l.sig <- length(unique(df.sig$Signatures))

df.sig$Cell.lines <- rep(colnames(cse), each =l.sig)
df.sig$Cell.lines <- factor(df.sig$Cell.lines, levels = unique(df.sig$Cell.lines)) 




# svg(filename = file.path(image.dir, paste0("de-novosiganturesCellLinesTest_",num.sign,".svg")), width = 13, height = 9)

ggplot(df.sig, aes(Cell.lines, Mutation.counts.by.signature, fill = Signatures)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, size = 20, hjust = 1, vjust = 0),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        panel.background = element_rect(fill = "white", color="white", size = 2),
        panel.grid= element_line(colour = "grey", size = 0.5, linetype = "solid") )
# dev.off()



msig1 <- matchSignatures(mutSign = cell.signature.sig,
                         reference = cosmic,
                         threshold = 0.45,
                         plot = T)
# png(file.path(image.dir,paste0("de-novoSigantures_MatchCosmic",num.sign, ".png")),width = 1000,height = 500)
msig1
# dev.off()


# After obtaining the match signatures results,
# # we select those cosmic signatures with sum coefficient Beta >0.2.
cell.sigs.cosmic <- mutSignatures::resolveMutSignatures(mutCountData = cells.signature.counts,
                                                        signFreqData=cosmic)

cell.sigs.cosmic.counts <- mutSignatures::as.data.frame(cell.sigs.cosmic$Results$count.result)
cscc <- t(cell.sigs.cosmic.counts)
cscc <- cscc[sample.names,]
cscc <- cell.sigs.cosmic.counts
write.table(cscc, file.path(mut.dir, "MutSignatures.SBS_AllSamples_ValSet.csv"), sep = "\t", col.names = TRUE, row.names = TRUE)

rownames(cscc)
# coef.beta <- mutSignatures::as.data.frame(cell.sigs.cosmic$results$freq.result)
# cosmic.selected<- rownames(coef.beta)[rowSums(coef.beta)>0.2]
# coef.beta <- coef.beta[cosmic.selected,]
# 
# png(file.path(image.dir,"Pheatmap_0.2CosmicSigantures_CoeffBeta.png"),width = 1000,height = 500)
# pheatmap(coef.beta)
# dev.off()
# # 
# cosmic.selected <- c("COSMIC.2","COSMIC.5","COSMIC.7","COSMIC.30")
# # cosmic.selected <- rownames(msig1$distanceMatrix)[rowSums(msig1$distanceMatrix)<0.6]
# cosmic.df <- mutSignatures::as.data.frame(cosmic)
# cosmic.df <- cosmic.df[,cosmic.selected]
# cosmic.selected <- mutSignatures::as.mutation.signatures(cosmic.df)
# # cosmic.selected <- cosmic

cell.sigs.cosmic <- mutSignatures::resolveMutSignatures(mutCountData = cells.signature.counts,
                                                        signFreqData=cosmic)

counts.cell.sig <- cell.sigs.cosmic$Results$count.result

png(file.path(image.dir,"SelectedCosmicsigantures_CellLines.png"),width = 1000,height = 500)
msigPlot(counts.cell.sig, sample =colnames(counts.cell.sig) )+ 
  scale_fill_manual(values = c(col_vector[1:80]))
dev.off()


# coef.beta <- mutSignatures::as.data.frame(cell.sigs.cosmic$results$freq.result)
# pheatmap(coef.beta)
# 
# signatures <- as.numeric(gsub("COSMIC.","", colnames(as.data.frame(cosmic.selected))))
# for(i in seq_len(length(signatures))){
#  
#   png(file.path(image.dir,paste0("COSMIC_", signatures[i],"_SelectedCosmicsigantures.png")),width = 1000,height = 500)
#   
#   msigPlot(cosmic.selected, signature = i, ylim = c(0, 0.3))
#   dev.off()
# }
# 
# 
# msig2 <- matchSignatures(mutSign = cell.signature.sig,
#                          reference = cosmic.selected, 
#                          threshold = 0.45,
#                          plot = T)
# png(file.path(image.dir,"de-novoSigantures_MatchSelectedCosmic.png"),width = 1000,height = 500)
# msig2
# dev.off()



cse <- mutSignatures::as.data.frame(counts.cell.sig)
#Filtering those signatures without counts
cse <- cse[-which(rowSums(cse)==0),]
colnames(cse)
for(i in seq_len(length(sample.names))){
  colnames(cse)[colnames(cse) %in% sample.names[i]] <- names(sample.names)[i] 
} 
pnf <- colnames(cse)[c(9,11,15,17,19,33)]
names(pnf) <- paste0("Blood-",seq(length(pnf)))
c1 <- c1[length(c1):1]
c2 <- c2[length(c2):1]
selected.samples <- c(pnf,c1,c2)
selected.samples <- c(paste0("BLOOD_",sample.data$Sample.name[sample.data$BAM_Normal]),sample.names)
cse <- cse [,selected.samples]
for(i in seq_len(length(selected.samples))){
  colnames(cse)[colnames(cse) %in% selected.samples[i]] <- names(selected.samples)[i] 
} 
#tumor wsg
if(stk.var[i]== "germ.variants"){
  # cse <- cse[,c(8,6,10,13,15,17,19,21,25,28,1,23,2,3,4,5,7,9,11,12,14,16,18,20,22,24,26,27,29)] 
  # cse <- cse[,c(7,9,12,14,16,26,1,2,6,3,4,5,34,11,13,15,17,18,19,21,23,24,25,27,30,32,33,35,20,8,10,22,29,28)]
  cse <- cse[,c(1,2,3,4,34,6,5,17,18,19:21,23,24,25,27,30,32,33,35,36,31,22,29,28)] #IB selected samples
  cse <- cse[,c()]
}
order.cse <- c(colnames(cse)[colnames(cse) %in%c1],colnames(cse)[colnames(cse) %in%c2])
names(order.cse) <- order.cse

order.cse <- order.cse[(sample.names)]
names(order.cse) <- sample.names
unique(colnames(cse)) %in% names(sample.names)


#tumor wgs
if(stk.var[i]== "germ.variants"){
  cse <- cse[,c(5,15,18,21,23,25,27,28,31,32,1:3,6:9,4,16,19,20,22,24,29,30,35,36,33,34,12:14,10,11)] #NoMPNST35
}

cse <- cse[,c(14,6,27,9,8)] #EdgarSignatures

df.sig <- c()
for(i in seq_len(ncol(cse))){
  sn <- colnames(cse)[i]
  df.sig<- c(df.sig,cse[,sn])
  
}


df.sig <- data.frame("Mutation counts by signature" = df.sig, COSMIC.Signatures= rep(rownames(cse),ncol(cse)))
df.sig$COSMIC.Signatures <- gsub("COSMIC.","",df.sig$COSMIC.Signatures)
l.cosmic <- length(unique(df.sig$COSMIC.Signatures))
df.sig$Cell.lines <- rep(colnames(cse), each =l.cosmic)
df.sig$Cell.lines <- factor(df.sig$Cell.lines, levels = unique(df.sig$Cell.lines)) 
df.sig$COSMIC.Signatures <- factor(df.sig$COSMIC.Signatures,levels= unique(df.sig$COSMIC.Signatures))

colnames(df.sig)
# df.sig$Mutation.counts.by.signature
set.seed(1214)
test.col <- sample(col_vector,nrow(cse))
test.col[1] <- "#6A9972"
test.col[6] <-"#8A0149"
test.col[15] <- "#595959"
test.col[17] <-"#943CF3"
i=2
# svg(filename = file.path(image.dir, paste0("de-novosigantures_MatchSelectedCosmic_",num.sign,".svg")), width = 13, height = 9)
# svg(filename = file.path(image.dir, paste0("Fitting_MatchAllCosmic_CellLines.v4.2.svg")), width = 13, height = 9)
svg(file.path("./MutSignatures/images", paste0(stk.var[i], "_Fitting_MatchAllCosmic_IntBiology_SelSamplesNames.v7.svg")), width = 25, height = 9)
# svg(file.path("./MutSignatures/images", paste0(stk.var[i], "_Fitting_MatchAllCosmic_Edgar_SelSamples.v6.svg")), width = 25, height = 15)

ggplot(df.sig, aes(Cell.lines, Mutation.counts.by.signature, fill = COSMIC.Signatures)) +
  geom_bar(stat="identity", position = "stack") + 
  # scale_fill_brewer(palette = "primary.colors")+
  scale_fill_manual(values = test.col)+
  theme(axis.text.x = element_text(angle = 90, size = 20, hjust = 1, vjust = 0, color = "black"),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 20,colour = "black"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        panel.background = element_rect(fill = "white", color="white", size = 2),
        panel.grid= element_line(colour = "grey", size = 0.5, linetype = "solid"))
dev.off()

order.nf.sp <- unique(df.sig$Cell.lines)
order.nf.sp <- order.nf.sp[c(1:3,8,4:7,14:27,9:13)]

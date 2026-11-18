#########################################
#       Preparation of Circos files     #
#########################################
library(regioneR)
######################## Functions ###########################

#If an error occurs, clean up as much as possible and exit with an error

changeTemplate <- function(x, template, replacement) {
  #message(paste0("Changing template   ", template, "   with   ", replacement ))
  x <- gsub(x=x, pattern=template, replacement=replacement, fixed=TRUE)
  return(x)
}
# First we need to generate a file per sample containing all ICR and INV
#INV and ICR will have an specific color.
execution.dir <-"."

sample.data <- read.table(file = file.path(execution.dir,"Sample.info.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
sample.names <- sample.data$Bioinfo.name[-c(1,30:32)]
names(sample.names) <- sample.data$Sample.name[-c(1,30:32)]

#Creating circos directory
for(i in seq_len(length(sample.names))){
  circos.dir <- file.path(execution.dir,"results",names(sample.names)[i],"Circos")
  if(!file.exists(circos.dir))dir.create(circos.dir)
}

#Load your procesed data
# load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.v5.RData")))



circos.frst.data <-list()  
# i=1
icr.col <- "color=(83,70,219)"
inv.col <- "color=(219,142,69)"
# filt.list$GFS_MPNST$AllRegions
for(i in seq_len(length(filt.list))){
  sn <- names(filt.list)[i]
  circos.data <- filt.list[[sn]]$AllRegions
  if(length(circos.data)==0)next
  
  circos.data <- circos.data[circos.data$SVTYPE %in%c("ICR","INV")]
  if(length(circos.data)==0)next
  mcols(circos.data) <-mcols(circos.data)[c(6,4,13)]
  # cliff <- join.cliff.lumpy[[sn]]
  # cliff <- cliff[cliff$SVTYPE%in%c("ICR","INV")]
  # mcols(cliff) <- mcols(cliff)[c(1,2,3)]
  # circos.data <- sort(unique(c(lumpy.icr.inv,cliff)))
  
  circos.data$SVTYPE[circos.data$SVTYPE == "ICR"] <- icr.col
  circos.data$SVTYPE[circos.data$SVTYPE == "INV"] <- inv.col
  circos.data.start <- data.frame(chr=as.character(seqnames(circos.data)), start= start(circos.data),end =end(circos.data))
  circos.data.start$chr <- gsub("chr","hs",circos.data.start$chr)
  circos.data.end <- toDataframe(toGRanges(circos.data$END))
  circos.data.end[,1] <- gsub("chr","hs",circos.data.end[,1])
  circos.data <- cbind(circos.data.start, circos.data.end,circos.data$SVTYPE)
  if(!file.exists( file.path(execution.dir,"Circos_files")))dir.create(file.path(execution.dir,"Circos_files"))
  # write.table(circos.data,file = file.path(execution.dir,"Circos_files",paste0("translocations_",sn,".txt")), sep="\t", col.names = FALSE,row.names = FALSE, quote = FALSE)
  # write.table(circos.data,file = file.path(execution.dir,"Circos_files",paste0("translocations_",sn,".txt")), sep="\t", col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(circos.data,file = file.path(execution.dir,"Circos_files",paste0("translocations_nofilt_",sn,".txt")), sep="\t", col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  
}

########## Modification of circos.transloc.conf #######
test <- "circos.transloc.nofilt"
# test <- "circos.transloc"

# sample.names <- sample.names$Sample.name[-c(1,5,30,31,32)]
for (i in seq_len(length(sample.names))){
  sample.name <-names(sample.names)[i]
  # cluster.file <- file.path(execution.dir, "toRunInCluster", paste0("Run_in_cluster_",test, ".sh"))
  cluster.file <- file.path(execution.dir,"Circos_files", paste0(test, ".conf"))
  # cluster.file <-'./Circos_files/circos.transloc.conf'
  # t <- thr[sample.name,]
  # fq.dir <- sample.data$file.fastq[i]
  content <- readLines(cluster.file)
  # content <- changeTemplate(content, "{%SAMPLE_DIR}", execution.dir)
  content <- changeTemplate(content,"{%SAMPLE_NAME%}", sample.name)
  new.cluster.file <- file.path("./Circos_files/", paste0(test,"_", sample.name,"_Nofilt.conf"))
  
  write(content, file = new.cluster.file)
  print(paste0("runcircos --conf=", new.cluster.file))
  ### Lunch pipeline
  # if (test == "StrelkaGerm"){
  #   # if (length( list.files(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_Tumors/WGS/results/", sample.name,"/Strelka_GermVariants"))) != 0) next
  # }else{
  #   if(length(list.files(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_Tumors/WGS/results/", sample.name,"/Manta_ManConfig_hg38"))) != 0) next
  # }
  
  system(command = paste0("runcircos --conf=", new.cluster.file),wait = TRUE )
  
}

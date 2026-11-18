#######################################################################
#                        Salmon Alignment                             #
#######################################################################
library(org.Hs.eg.db)
####### Parameters
analysis.dir <- "."
sample.data <- read.table(file = "./Sample.info.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE,comment.char = "")
file.names <- sample.data$Sample.name

# Salmon alignement and quantification parameters
# salmonDir <- "salmon"
salmonDir <- "./salmon-1.8.0/bin/salmon"
file1.suffix <- "_1.fastq.gz"
file2.suffix <- "_2.fastq.gz"
fastqdir <- "./Data/"
transcript.index <- "./Salmon_index/salmon_indexes_UCSC_hg38"
output.suffix <- "_quant"
output.quants <- "./results/Salmon/"
threads <- 6
source(file = file.path("./RNA-seq_MesennchimalDiff/utils.R"))

### salmonAlignment
salmonAlignment <- function(sample.name, salmonDir,
                            file1.suffix, file2.suffix,
                            fastqdir,
                            transcript.index,
                            output.suffix,
                            output.quants, 
                            libtype = "IU",
                            threads = 4,
                            verbose = TRUE){
  now.msg("Salmon starting...")
  file.name <- file.path(output.quants, paste0(sample.name, output.suffix))
  
  if (!file.exists(file.name)){
    now.msg("missing file", file.name, ". Obtaining quant file")
    full.command <- paste0(salmonDir, " quant -i",
                           transcript.index," -l ",
                           libtype, " -1 ", 
                           fastqdir,
                           sample.name,
                           file1.suffix, " -2 ", 
                           fastqdir, sample.name, 
                           file2.suffix, " --validateMappings -p ",
                           as.character(threads),
                           " -o ", 
                           output.quants,
                           sample.name, output.suffix )
    system(full.command, wait = TRUE)
  }else{
    now.msg("Salmon quants already obtained")
  }
}


###################          Salmon alignment           #####################

# ## Executing Salmon alignement by Selective allignement and quantification
# sn <- names(file.names)[1]
file.names <- "54MPNSTB"
for(i in seq_len(length(file.names))){
  sn <- file.names[i]
  salmonAlignment(sample.name = sn, salmonDir,
                  file1.suffix = file1.suffix,
                  file2.suffix = file2.suffix,
                  fastqdir = fastqdir,
                  transcript.index = transcript.index,
                  output.suffix = output.suffix,
                  output.quants = output.quants,
                  threads = threads)
}
now.msg("   Salmon done")

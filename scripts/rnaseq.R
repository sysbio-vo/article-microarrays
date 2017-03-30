library(Rsubread)
library(limma)

targets <- readTargets()

INDEX = TRUE

if (INDEX) {
  buildindex(basename="../hg19/hg19_index",reference="../hg19/hg19.fa")
  INDEX=FALSE
}

align(index="../hg19/hg19_index", unique=TRUE, indels=5,
      input_format="gzFASTQ", output_format="BAM", nthreads = 20, 
      readfile1=paste("../raws/", targets$InputFile1[1], sep=""),
      readfile2=paste("../raws/", targets$InputFile2[1], sep=""),
      output_file=paste("../output/", targets$OutputFile[1], sep=""))

#featureCounts
#voom


library(Rsubread)
library(limma)

targets <- readTargets()

GSE = "GSE58135"
INDEX = TRUE

if (INDEX) {
  buildindex(basename="../hg19/hg19_index",reference="../hg19/hg19.fa")
  INDEX=FALSE
}

align(index="../hg19/hg19_index", unique=TRUE, indels=5,
      input_format="gzFASTQ", output_format="BAM", nthreads = 20, 
      readfile1=paste("../raws/", GSE, "/", targets$InputFile1[2:10], sep=""),
      readfile2=paste("../raws/", GSE, "/", targets$InputFile2[2:10], sep=""),
      output_file=paste("../bam/", GSE, "/", targets$OutputFile[2:10], sep=""))

#featureCounts
#voom


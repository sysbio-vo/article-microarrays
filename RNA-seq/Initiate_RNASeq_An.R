####################################################################################################################################
############################################ Initiate RNA-Seq data analysis ########################################################
####################################################################################################################################
#Pass arguments:
args<-commandArgs(trailingOnly = TRUE)
start<-as.numeric(as.character(args[1]))
end<-as.numeric(as.character(args[2]))
# Read pheno data
data_pheno<-read.delim("pdata_GSE58135.txt",sep="\t",header=F)
data_pheno<-data_pheno[which(data_pheno$V12=="tissue: Triple Negative Breast Cancer Primary Tumor"),]
file_dirs<-as.character(data_pheno$V43)
#Define subset of samples:
file_dirs<-file_dirs[start:end]
for(file in file_dirs){
  sub<-unlist(strsplit(file,"/"))
  SRX<-sub[length(sub)]
  command<-paste("run.sh -v SRX=",SRX,sep="",collapse = "")
  system(command) 
}
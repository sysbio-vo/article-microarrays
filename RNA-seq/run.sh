#DIR="/home/user/article-microarrays/RNASeq/"
SRA=".sra"
FQ_1="_1.fastq"
FQ_2="_2.fastq"
OUTt="_tophat"
OUTc="_cufflinks"

################ Run RNA-Seq data processing
#cd $DIR
# Download data:
wget -r --no-parent ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX559/$SRX
# Reformat:
SDIR="ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX559/"$SRX
SRR="$(ls $SDIR)"
mv $SDIR/$SRR/* ./
# Unpack with fastq-dump:
fastq-dump.2.5.0 --split-files $SRR$SRA -O ./
# Run TopHat:
tophat2 hg19 -p 5 -o $SRR$OUTt -G RefSeq_hg.gtf $SRR$FQ_1,$SRR$FQ_2
# Run Cufflinks:
cufflinks -p 5 -G RefSeq_hg.gtf -o $SRR$OUTc $SRR$OUTt/accepted_hits.bam













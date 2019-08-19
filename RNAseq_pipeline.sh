##### Pipeline Information for Alex Waugh 

#### This is a loose pipeline for the analysis of RNA-seq data from raw reads to differential expression. All codes are written to run on the sapelo2 cluster

### First, where's the data

/project/bghlab/Sinv_Overwinter
# only accessible through qlogin and xfer nodes. I recommend moving the .fq files to your scratch directory to work with
# Files beginning with "SF" are "Spring Flying" and "OW" are "Overwintered". 
# "Br" indicates brain tissue, "Ov" indicates "Ovary"
# Individual IDs are given by the number. OW_Br1.fq and OW_Ov1.fq are from the same individual (Though not the same as SF_Br1.fq and SF_Ov1.fq)

### Quality Control Steps

## Program: Fastqc

module load FastQC/0.11.8-Java-1.8.0_144

run_fastqc () {
 fastqc --extract -o ./fastqc_out $1
}

run_fastqc OW_Br1.fastq

## Program: Trim_Galore

module load Trim_Galore/0.4.5-foss-2016b

run_trimGalore () {
trim_galore --fastqc --gzip --paired ../simple_data/$1_R1.fq.gz ../simple_data/$1_R2.fq.gz
}

run_trimGalore G_pBL_1

## Check post-trimming fastqc reports. Might need to remove rRNA but probably unnecessary
# If needed, check with Sam and he can provide a script that's usable for this. 

### Alignment Steps

## Program: STAR
# Genome Preparation

ml STAR/2.7.1a-foss-2016b

STAR --runMode genomeGenerate --runThreadN 8 --genomeSAindexNbases 13 --genomeDir /scratch/sva/LK_analysis/genome/STAR_genome_GTF --genomeFastaFiles /scratch/sva/LK_analysis/genome/SINVBB1.nr.fa --sjdbGTFfile /scratch/sva/LK_analysis/genome/SINVBB1.chr.gtf

# Note: Annotation needs to be converted to gtf if necessary
module load Cufflinks/2.2.1-foss-2016b
gffread my.gff3 -T -o my.gtf

# Also, pay attention to error messages. May need to adjust --genomeSAindexNbases depending on genome. 

## Program STAR
# Alignment first-pass
run_STAR_zipped () {
time STAR --readFilesCommand zcat --runThreadN 4 --genomeDir /scratch/sva/LK_analysis/genome/chr16_genome/ --readFilesIn /scratch/sva/LK_analysis/trimmed_reads/$1_trimmed.fq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFileNamePrefix /scratch/sva/LK_analysis/smallb_aln/$1. --limitBAMsortRAM 15000000000 
}
run_STAR_unzipped () {
time STAR --runThreadN 4 --genomeDir /scratch/sva/LK_analysis/genome/chr16_genome/ --readFilesIn /scratch/sva/LK_analysis/gyne_data_trim/$1_trim.fq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFileNamePrefix /scratch/sva/LK_analysis/smallb_aln/$1. --limitBAMsortRAM 15000000000 
}

run_STAR_unzipped 104GB

## Align libraries again with the addition of the junction information passed to the aligner using the --sjdbFileChrStartEnd parameter and the --quantMode GeneCounts parameter to generate gene counts for each libary

### Differential Expression Steps
## Program: EdgeR

# Download the counts files to a directory on your laptop and install the edgeR package in R. I generally recommend using R through RStudio which is a free IDE that is fantastic for scripting. 
## Here is an example of a differential expression analysis code that I wrote for reference. I recommend reading through the edgeR documentation


## Clear Workspace
rm(list=ls())

## Import Libraries
library(edgeR) 
## Load in functions

## Load Files individually into a data.frame 
## Alternatively, you can merge the gene counts files manually or using shell scripting if you prefer that method and load them in to R directly. 
setwd("~/Desktop/LK_DiffExp/Worker_DE")
temp = list.files(pattern="*.ReadsPerGene.out.tab")
temp2 <- sapply(strsplit(temp, split='.', fixed=TRUE), function(x) (x[1]))
myfiles = lapply(temp, read.delim,header=F)
myfiles_new <- myfiles
# Note: Each File is saved as myfiles[[x]] as a data.frame

## merge first two columns into identifier column
for (i in 1:26){
  if (grepl("B_",temp2[i])) {
    colnames(myfiles_new[[i]])[2] <- paste(temp2[i],"_counts",sep="")
    myfiles_new[[i]] <- myfiles_new[[i]][,c(1,2)]
  } else {
    colnames(myfiles_new[[i]])[4] <- paste(temp2[i],"_counts",sep="")
    myfiles_new[[i]] <- myfiles_new[[i]][,c(1,4)]
  }
}

## merge all data.frames by ID and reformat into gene_data.frame
merged <- Reduce(function(...) merge(..., by='V1', all.x=FALSE), myfiles_new)
row.names(merged) <- merged$V1
x <- merged[,2:27]

## Load in Model Matrix
adv_group <- read.delim("/Users/samarsenault/Desktop/LK_DiffExp/Worker_DE/Sinv_model_matrix_worker.tsv",row.names = 1,sep = "\t")
## Model Matrices generally look something like this:
"
Sample	Individual_ID	Colony_ID	Tissue	Social Form	Genotype	Alt_ID
104AB	104A	104	Brain	Polygyne	Bb	Brain.Polygyne.Bb
104AO	104A	104	Ovary	Polygyne	Bb	Ovary.Polygyne.Bb
104DB	104D	104	Brain	Polygyne	BB	Brain.Polygyne.BB
104GB	104G	104	Brain	Polygyne	bb	Brain.Polygyne.bb
104GO	104G	104	Ovary	Polygyne	bb	Ovary.Polygyne.bb
"
## They provide relevant information on the nature of the samples so that the generalized linear models edgeR uses can account for that variation when searching for differences due to specific factors.
 
## Pairwise Comparisons
y <- DGEList(x,group = adv_group$Alt_ID)
design_pair <- model.matrix(~0+adv_group$Alt_ID,data=y$samples)
keep <- rowSums(cpm(y)>1) >= 6 ## Troubleshoot this parameter. It eliminates low coverage genes that you won't be able to get high-confidence results for. 
y <- y[keep, , keep.lib.sizes=FALSE]
colnames(design_pair) <- levels(y$samples$group)

## Normalize samples by library size
y <- calcNormFactors(y)
Background.genes <- row.names(y$counts)
y <- estimateDisp(y,design_pair)
y <- estimateGLMCommonDisp(y, design_pair)
y <- estimateGLMTrendedDisp(y, design_pair)
y <- estimateGLMTagwiseDisp(y, design_pair)
fit <- glmQLFit(y,design_pair)

my.contrasts <- makeContrasts(Brain_BB_v_Bb=Brain_pBB-Brain_pBL, 
                              Gaster_BB_v_Bb=Gaster_pBB-Gaster_pBL,
                              levels=design_pair)

## Genotypic Changes

qlf.Brain_BB_v_Bb <- topTags(glmQLFTest(fit, contrast=my.contrasts[,"Brain_BB_v_Bb"]),p.value = 1,n=100000)$table
qlf.Gaster_BB_v_Bb <- topTags(glmQLFTest(fit, contrast=my.contrasts[,"Gaster_BB_v_Bb"]),p.value = 1,n=100000)$table


save.image("~/Desktop/LK_DiffExp/Worker_DE/Output.RData") ## This helps alot if the code you're running takes a while and you don't want to have to keep running it. Save your workspace and you can restart from there. 
load("~/Desktop/LK_DiffExp/Worker_DE/Output.RData")

y_cpm <- cpm(y)
save(y_cpm,file ="~/Desktop/LK_DiffExp/Worker_DE/Worker_CPM.RData" ) ## Saves a particular R variable to be loaded into other R programs. 
write.table(y_cpm,file =  "~/Desktop/LK_DiffExp/Worker_DE/Worker_CPM.tsv",sep = "\t",row.names=TRUE,col.names = TRUE,quote = FALSE) ## Writes the results into a more generally readable format. 





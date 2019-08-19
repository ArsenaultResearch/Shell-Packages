##### Bisulfite Sequencing Analysis 2019

#### This is a loose pipeline for the analysis of BS-seq data from Raw reads on the sapelo2 cluster

### Description of the data

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

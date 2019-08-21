##### Bisulfite Sequencing Analysis 2019

#### This is a loose pipeline for the analysis of BS-seq data from Raw reads on the sapelo2 cluster

### Description of the data


### Quality Control Steps


## Program: Trim_Galore

module load Trim_Galore/0.4.5-foss-2016b

run_trimGalore () {
trim_galore --fastqc --gzip --paired ../simple_data/$1_R1.fq.gz ../simple_data/$1_R2.fq.gz
}

run_trimGalore G_pBL_1

## Check post-trimming fastqc reports. Might need to remove rRNA but probably unnecessary
# If needed, check with Sam and he can provide a script that's usable for this. 


### Prepare the genome for alignment

## Program: Bismark
module load Bismark/0.22.1-foss-2016b
module load Bowtie2/2.3.4.1-foss-2016b

bismark_genome_preparation --verbose /scratch/sva/Sinv_Methylation/genome/


### Align the reads to the genome

## Program Bismark
module load Bismark/0.22.1-foss-2016b
module load Bowtie2/2.3.4.1-foss-2016b
module load SAMtools/1.9-foss-2016b


run_bismark () { 
bismark --genome /scratch/sva/Sinv_Methylation/genome/ --fastq --gzip --sam --output_dir /scratch/sva/Sinv_Methylation/analyses/bismark_aln --temp_dir /scratch/sva/Sinv_Methylation/analyses/bismark_temp --basename $1 -1 /scratch/sva/Sinv_Methylation/analyses/trimmed_reads/$1_R1_val_1.fq.gz -2 /scratch/sva/Sinv_Methylation/analyses/trimmed_reads/$1_R2_val_2.fq.gz
deduplicate_bismark $1_pe.sam
bismark_methylation_extractor 
}

### Generate Reports from the bismark outputs
## Program: Bismark
run_bismark_reporting () { 
bismark2report
bismark2summary
}
run_bismark_reporting

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
bismark --multicore 3 --genome /scratch/sva/Sinv_Methylation/genome/ --fastq --gzip --output_dir /scratch/sva/Sinv_Methylation/analyses/bismark_aln --temp_dir /scratch/sva/Sinv_Methylation/analyses/bismark_temp -1 /scratch/sva/Sinv_Methylation/analyses/trimmed_reads/$1_R1_val_1.fq.gz -2 /scratch/sva/Sinv_Methylation/analyses/trimmed_reads/$1_R2_val_2.fq.gz
deduplicate_bismark --bam -p $1_R1_val_1_bismark_bt2_pe.bam
bismark_methylation_extractor --paired-end --CX --cytosine_report --genome_folder /scratch/sva/Sinv_Methylation/genome/ $1_R1_val_1_bismark_bt2_pe.deduplicated.bam
}
## Note: run this on a 9-core run. Bismark mentioned that they require 3 cores for every 1 value in multicore. 
## --multicore is incompatible with --basename at the moment which is slightly inconvenient but may change in the near future. 

### Generate Reports from the bismark outputs
## Program: Bismark
run_bismark_reporting () { 
bismark2report
bismark2summary
}
run_bismark_reporting

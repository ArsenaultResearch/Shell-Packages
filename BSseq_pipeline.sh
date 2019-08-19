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

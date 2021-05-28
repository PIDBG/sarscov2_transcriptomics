- Ali's iBSc project: mining publicly available SARS-CoV-2 transcriptomic datasets to analyse severity and validate signatures

##### Build index genome (GRCH38 assembly) using STAR
#PBS -l select=1:ncpus=24:mem=60gb
#PBS -l walltime=24:00:00


module load anaconda3/personal 
source activate /rds/general/user/ah2720/home/anaconda3

STAR --runThreadN 23 --runMode genomeGenerate --genomeDir /rds/general/user/ah2720/home/STAR3 --genomeFastaFiles /rds/general/user/ah2720/home/genome/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /rds/general/user/ah2720/home/genome/gencode.v37.annotation.gtf


##### Download FASTQ files

#PBS -l select=1:ncpus=32:mem=62gb
#PBS -l walltime=24:00:00


export PATH=$PATH:/rds/general/user/ah2720/home/sradownloader-3.6
export PATH=$PATH:/rds/general/user/ah2720/home/sratoolkit.2.10.9-ubuntu64/bin

cd /rds/general/user/ah2720/projects/covid_diamonds/live/GSE_161731_new
sradownloader --outdir fastq_161731 SraAccList_161731_1.txt


##### Trim FASTQ files - trimmomatic

#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=72:00:00

module load trimmomatic

cat /rds/general/project/covid_diamonds/live/GSE_161731_new/SraAccList_161731_QC_part1.txt | while read id
do
trimmomatic PE -threads 8 -phred33 \
 /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731/$id*_1.fastq.gz /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731/$id*_2.fastq.gz \
 /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731_trimmed/$id*_1.fastq.gz /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731_trimmed_rubbish/$id*_unpaired_1.fastq.gz \
 /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731_trimmed/$id*_2.fastq.gz /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731_trimmed_rubbish/$id*_unpaired_2.fastq.gz \
 ILLUMINACLIP:/rds/general/project/covid_diamonds/live/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


##### Map to index Genoma (pre-built)

#!/bin/bash
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=48:00:00

module load anaconda3
source activate /rds/general/user/ah2720/home/anaconda3

cat /rds/general/project/covid_diamonds/live/GSE_161731_new/SraAccList_161731_QC_part1.txt | while read id
do
STAR \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--readFilesCommand zcat \
--runThreadN 8 \
--sjdbGTFfile /rds/general/user/ah2720/home/genome/gencode.v37.annotation.gtf \
--outReadsUnmapped Fastx \
--outMultimapperOrder Random \
--outWigType wiggle \
--genomeDir /rds/general/user/ah2720/home/STAR3 \
--readFilesIn /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731_trimmed/$id*_1.fastq.gz /rds/general/project/covid_diamonds/live/GSE_161731_new/fastq_161731_trimmed/$id*_2.fastq.gz \
--outFileNamePrefix /rds/general/project/covid_diamonds/live/GSE_161731_new/aligned_161731_new/$id 
done


##### Count BAM files using FeatureCounts

#PBS -l select=1:ncpus=32:mem=62gb
#PBS -l walltime=72:00:00

module load anaconda3
source activate /rds/general/user/ah2720/home/anaconda3

featureCounts -T 30 -p -s 1 -a /rds/general/user/ah2720/home/genome/gencode.v37.annotation.gtf -g gene_id -o /rds/general/project/covid_diamonds/live/GSE_161731_new/readCounts_new_software_update/readCounts_new.txt /rds/general/project/covid_diamonds/live/GSE_161731_new/aligned_161731_new/*.bam

##### Cleaning up FeatureCounts matrix - Linux

https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html --> link used


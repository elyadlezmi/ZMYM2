#!/usr/bin/bash

## RNA-seq pipeline
p=8

single="emp1 emp2 emp3 z6 z15 z22"

for file in $single;

do

/home/elyad/snap/TrimGalore-0.6.0/trim_galore -J $p $file.fastq

STAR --runThreadN $p --genomeDir /home/elyad/reference_genome/STAR --readFilesIn ${file}_trimmed.fq --quantMode GeneCounts --outFileNamePrefix ./$file --outSAMtype BAM SortedByCoordinate

rm -rf $file.fastq_trimming_report.txt ${file}_trimmed.fq ${file}Log.out ${file}Log.progress.out ${file}SJ.out.tab 

samtools index ${file}Aligned.sortedByCoord.out.bam

done



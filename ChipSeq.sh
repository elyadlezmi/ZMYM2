#!/usr/bin/bash

## chip-seq pipeline

#rename 's/_S\w+001//' *

p=8

for file in emp1_H3Ac emp2_H3Ac emp3_H3Ac z6_H3Ac z15_H3Ac z22_H3Ac emp1_H3K4me1 emp2_H3K4me1 emp3_H3K4me1 z6_H3K4me1 z15_H3K4me1 z22_H3K4me1;

do

echo $file

#bowtie2 -p $p -x /home/elyad/LabSync/ChipReference/GRCh38 -U ${file}.fastq.gz | samtools view -b -q 5 | samtools sort > $file.undup.bam

#samtools rmdup -s $file.undup.bam $file.bam
#rm -rf $file.undup.bam

#samtools index $file.bam

#bamCoverage -b $file.bam --outFileFormat bigwig -o ${file}.bw -p $p --normalizeUsing CPM --smoothLength 150 --exactScaling --effectiveGenomeSize 2913022398

#bamCoverage -b $file.bam --outFileFormat bedgraph -o ${file}.bed -p $p --normalizeUsing CPM --exactScaling --effectiveGenomeSize 2913022398

macs2 callpeak -B -t $file.bam -n $file --keep-dup 1 -g hs --broad --nomodel --extsize 400 --broad-cutoff 0.05

#macs2 callpeak -B -t $file.bam -n $file --keep-dup 1 -g hs -q 0.1

rm -rf ${file}_peaks.gappedPeak ${file}_treat_pileup.bdg ${file}_peaks.xls ${file}_peaks.gappedPeak ${file}_summits.bed ${file}_control_lambda.bdg ${file}_model.r

done



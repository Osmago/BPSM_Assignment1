#!bin/bash

#this script will get data from RNA sequencing of two life stages of Tripanosoma bruscei bruscei,
#with three samples each, and count the mean count of reads aligning with each reference gene

#first set the variables relating to each initial file

#fqfiles will guide a loop later to the paired reads files
guide="/localdisk/data/BPSM/Assignment1/fastq/fqfiles"

#genome file must be decompressed to be used by bowtie2
zip_genome="/localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz"
zcat $zip_genome > Tb927_genome.fasta
genome="Tbb_genome.fasta"

#genes file to generate counts in the end
genes="/localdisk/data/BPSM/Assignment1/Tbbgenes.bed"

#make or confirm location for fastqc outputs, then clean it
if test ! -d fastqc_outputs; then
 mkdir fastqc_outputs
 fi
#rm -fr fastqc_outputs/*

declare what are the files being used
echo -e "\nUsing $guide\tfile as guide for paired reads"
echo -e "Using $zip_genome\tfile as genome"
echo -e "Using $genes\tfile as reference for gene locations\n"

#sorts fqfiles just so we do everything in order, then pipes into data processing loop
sort -k1,1n $guide | while read smpl_nmbr smpl_type pair1 pair2; do
 pair1_path="/localdisk/data/BPSM/Assignment1/fastq/$pair1"
 pair2_path="/localdisk/data/BPSM/Assignment1/fastq/$pair2"
# echo -e "\nStarting read quality assessment of sample $smpl_nmbr\n"
# fastqc --extract -o fastqc_outputs/ $pair1_path $pair2_path
 echo -e "Results for quality assessmentof pair $smpl_nmbr:"
 echo -e "\n$pair1"
 cut -f1,2 fastqc_outputs/"$smpl_nmbr"_L8_1_fastqc/summary.txt
 echo -e "\n$pair2"
 cut -f1,2 fastqc_outputs/"$smpl_nmbr"_L8_2_fastqc/summary.txt
done

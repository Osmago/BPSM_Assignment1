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

#declare what are the files being used
echo "using $guide file as guide for paired reads"
echo "using $zip_genome file as genome"
echo "using $genes file as reference for gene locations"

while read smpl_nmbr smpl_type pair1 pair2; do
 pair1_path="/localdisk/data/BPSM/Assignment1/fastq/$pair1"
 pair2_path="/localdisk/data/BPSM/Assignment1/fastq/$pair2"
 
done < $guide

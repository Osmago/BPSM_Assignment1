#!bin/bash

#this script will get data from RNA sequencing of two life stages of Tripanosoma bruscei bruscei,
#with three samples each, and count the mean count of reads aligning with each reference gene

#first set the variables relating to each initial file

#fqfiles will guide a loop later to the paired reads files
guide="/localdisk/data/BPSM/Assignment1/fastq/fqfiles"

#genome file must be decompressed to be used by bowtie2
zip_genome="/localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz"
zcat $zip_genome > Tb927_genome.fasta
genome="Tb927_genome.fasta"

#genes file to generate counts in the end
genes="/localdisk/data/BPSM/Assignment1/Tbbgenes.bed"

#make or confirm location for fastqc outputs, then clean it
if test ! -d fastqc_outputs; then
 mkdir fastqc_outputs
 fi
rm -fr fastqc_outputs/*

#declare what are the files being used
echo -e "\nUsing $guide\tfile as guide for paired reads"
echo -e "Using $zip_genome\tfile as genome"
echo -e "Using $genes\tfile as reference for gene locations\n"

#create bowtie2 database for our genome as it will be used for every read
#omitted bowtie2 output because it's too big and the user should know the quality of the genome
genome_prefix=$(echo "$genome" | cut -f1 -d '.')
echo -e "Creating $genome_prefix bowtie2 database...\n"
bowtie2-build $genome $genome_prefix > /dev/null
echo -e "\nGenome database creation done"

#sorts fqfiles just so we do everything in order, then pipes into data processing loop
sort -k1,1n $guide | while read smpl_nmbr smpl_type pair1 pair2; do

 #saving pairs files names so i can use them later
 pair1_path="/localdisk/data/BPSM/Assignment1/fastq/$pair1"
 pair2_path="/localdisk/data/BPSM/Assignment1/fastq/$pair2"

 #fastqc analysis
 echo -e "\nAssessing sequencing quality of sample $smpl_nmbr...\n"
 fastqc --extract -o fastqc_outputs/ $pair1_path $pair2_path

 #fastqc results presented for user
 echo -e "\nResults for quality assessment of pair $smpl_nmbr:"
 echo -e "\n$pair1"
 cut -f1,2 fastqc_outputs/"$smpl_nmbr"_L8_1_fastqc/summary.txt
 echo -e "\n$pair2"
 cut -f1,2 fastqc_outputs/"$smpl_nmbr"_L8_2_fastqc/summary.txt
 echo
 
 #user can choose to use this sequencing data or not. if user doesn't type y, n, Y or N, it will ask again 
 chose=0
 while test $chose -eq 0; do
  read -p "Do you want to use this sequencing data? [y/n]" -n 1 yn </dev/tty
  case $yn in
   #if user chose yes, continue processing
   [Yy])
    chose=1

    #use bowtie2 pair-based sequencing options to align
    #using 3 threads to make it a little bit faster
    echo -e "\nPairing sample $smpl_nmbr reads to genome...\n"
    bowtie2 -x $genome_prefix -1 $pair1_path -2 $pair2_path -S $smpl_nmbr.sam --threads 3
    
    #use samtools to convert the output to BAM format ('view' gets SAM input and -b scpecifies outputs to BAM)
    samtools view -b $smpl_nmbr.sam > $smpl_nmbr.bam

    #use bedtools to get gene count for this read pair
    #pairtobed gets the hits for this read pair in this bed file
    #using -abam to specify BAM input and -bedpe to specify BED formatted output, where the gene name is in column 14. We will sort and count each gene hit
    echo -e "\nCounting read alignments per gene for sample $smpl_nmbr"
    bedtools pairtobed -abam $smpl_nmbr.bam -b $genes -bedpe | cut -f14 | sort | uniq -c > $smpl_nmbr.bed

    echo -e "\nNext read pair..."
   ;;
   #if user chose no, go to next pair 
   [Nn])
    chose=1
    echo -e "\nNext read pair..."
   ;;
   #if user chose anything else other than y,Y,n or N, ask again
   *)
    echo -e "\nPlease answer with y/n"
  esac
 done
done

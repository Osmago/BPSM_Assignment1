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

#define where to store which samples will be used by the user, in case she chooses not to use any sample. this will be useful for calculating means later
samples_used="samples_used"

#declare what are the files being used
echo
echo -e "Using $guide\tfile as guide for paired reads"
echo -e "Using $zip_genome\tfile as genome"
echo -e "Using $genes\tfile as reference for gene locations"
echo

#create bowtie2 database for our genome as it will be used for every read
#omitted bowtie2 output because it's too big and the user should know the quality of the genome
#using 10 threads to make it a little bit faster
genome_prefix=$(echo "$genome" | cut -f1 -d '.')
echo "Creating $genome_prefix bowtie2 database..."
echo
#checkpoint variable checks program run time before and after long processes, so that i can echo the processe's run time
checkpoint=$SECONDS
bowtie2-build $genome $genome_prefix --threads 10 > /dev/null
echo
echo "Genome database creation done"
echo "$((($SECONDS-$checkpoint) / 60)) minutes and $((($SECONDS-$checkpoint) % 60)) seconds elapsed."

#fastqc analysis
#using 10 threads to make it a little bit faster
#get a space-delimited list of the paths to the data of each 
fqpaths=$(awk '{a="/localdisk/data/BPSM/Assignment1/fastq/"; b=a$3" "a$4; print b}' $guide | paste -d ' ' -s)
echo
echo "Assessing sequencing quality of samples..."
echo
checkpoint=$SECONDS
fastqc --extract -o fastqc_outputs/ $fqpaths --threads 10
echo "Sequence quality assessment done"
echo "$((($SECONDS-$checkpoint) / 60)) minutes and $((($SECONDS-$checkpoint) % 60)) seconds elapsed."
echo

#this loop will separate the analysis of each sequencing pair, sample by sample
#this includes coosing wether to use the samples or not, based on the fastqc report summary
#then analysing the pair's data with bowtie2, samtools and bedtools
while read smpl_nmbr smpl_type pair1 pair2; do

 #saving pairs files names so i can use them later
 pair1_path="/localdisk/data/BPSM/Assignment1/fastq/$pair1"
 pair2_path="/localdisk/data/BPSM/Assignment1/fastq/$pair2"

 #fastqc results presented for user
 echo
 echo "Results for quality assessment of pair $smpl_nmbr:"
 echo
 echo "$pair1"
 cut -f1,2 fastqc_outputs/"$smpl_nmbr"_L8_1_fastqc/summary.txt
 echo
 echo "$pair2"
 cut -f1,2 fastqc_outputs/"$smpl_nmbr"_L8_2_fastqc/summary.txt
 echo
 
 #user can choose to use this sequencing data or not. if user doesn't type y, n, Y or N, it will ask again 
 chose=0
 while test $chose -eq 0; do
#  read -p "Do you want to use this sequencing data? [y/n]" -n 1 yn </dev/tty
  echo
  #auto-selecting yes to quicken testing process
  yn="y"
  case $yn in
   #if user chose yes, continue processing
   [Yy])
    chose=1

    #record that this is a sample being used
    echo -e "$smpl_nmbr\t$smpl_type" >> $samples_used

    #use bowtie2 pair-based sequencing options to align
    #using 10 threads to make it faster
    echo "Pairing sample $smpl_nmbr reads to genome..."
    echo
    checkpoint=$SECONDS
    bowtie2 -x $genome_prefix -1 $pair1_path -2 $pair2_path -S $smpl_nmbr.sam --threads 10
    echo "$((($SECONDS-$checkpoint) / 60)) minutes and $((($SECONDS-$checkpoint) % 60)) seconds elapsed."
    
    #use samtools to convert the output to BAM format ('view' gets SAM input and -b scpecifies outputs to BAM)
    #using 10 threads to make it a little bit faster
    echo
    echo "Preparing data for gene counts..."
    checkpoint=$SECONDS
    samtools view -@ 10 -b $smpl_nmbr.sam > $smpl_nmbr.bam

    #use bedtools to get gene count for this read pair
    #pairtobed gets the hits for this read pair in this bed file
    #using -abam to specify BAM input and -bedpe to specify BED formatted output, where the gene name is in column 14. We will sort and count each gene hit with sort and uniq
    echo
    echo "Counting read alignments per gene for sample $smpl_nmbr..."
    bedtools pairtobed -abam $smpl_nmbr.bam -b $genes -bedpe | cut -f14 | sort | uniq -c > $smpl_nmbr.counts
    echo
    echo "Finished counting sample $smpl_nmbr read alignments with genes"
    echo "$((($SECONDS-$checkpoint) / 60)) minutes and $((($SECONDS-$checkpoint) % 60)) seconds elapsed."
   ;;
   #if user chose no, go to next pair 
   [Nn])
    chose=1
    echo
    echo "Next read pair..."
   ;;
   #if user chose anything else other than y,Y,n or N, ask again
   *)
    echo
    echo "Please answer with y/n"
  esac
 done
done < <(sort -k1,1n $guide)

#now finally for the means calculation

#announce beggining of mean calculation
echo
echo "Calculating mean read counts for every gene among each of the sample types..."
checkpoint=$SECONDS

#define the output file (there are a lot of genes, we don't want to just show everything on the screen
final_output="Tbbgenes_means.tsv"

#clean previous outputs if there were any
rm -f $final_output

#do a header for the output
echo -e "gene_name\tSld_mean\tStp_mean" >> $final_output

#this loop will go through each gene in our genes file and check in the .counts files created before if there were counts for it
#if so, it will add each of the counts occurances and calculate a mean of all counts for each sample type in the end, appending this to the output file
while read gene; do

 #define variables used to count total number of counts and occurances in .counts files
 stumpy_total=0
 stumpy_counts=0
 slender_total=0
 slender_counts=0

 #this loop will go through each sample and add the value of the count number to a cumulative total, also counting how many samples were used
 while read smpl_nmbr smpl_type; do

  #the awk script will look for the current gene of the while loop in the .counts file and output its count number. outputs nothing if gene is not present
  counts=$(awk -v gene="$gene" '{if($NF == gene){print $(NF-1);}}' $smpl_nmbr.counts)

  #if statement prevents trying to add empty values
  if test "$counts" = ""; then
   counts=0
   fi

  #if statements add the count number to the cumulative total and count how many samples are being used, for each sample type
  if test "$smpl_type" = "Stumpy"; then
   stumpy_total=$(($stumpy_total+$counts))
   stumpy_counts=$(($stumpy_counts+1))
   fi
  if test "$smpl_type" = "Slender"; then
   slender_total=$(($slender_total+$counts))
   slender_counts=$(($slender_counts+1))
   fi

 done < $samples_used
 
 #if statements prevent dividing by zero if no Stumpy or Slender samples were used.
 #they then calculate the mean counts by dividing the cumulative total of read alignments by the number of samples used
 if test $stumpy_counts -ne 0; then
  stumpy_mean=$(($stumpy_total/$stumpy_counts))
 else
  stumpy_mean=0
  fi
 if test $slender_counts -ne 0; then
  slender_mean=$(($slender_total/$slender_counts))
 else
  slender_mean=0
  fi
 
 #finally append the mean values to our final output file
 echo -e "$gene\t$slender_mean\t$stumpy_mean" >> $final_output

#set the genes to be analysed by while loop by getting the 4th column of the genes file
done < <(cut -f4 $genes)

#cleanup samples used in this particular run of the program
rm -f $samples_used

#exit message
echo
echo "Mean values for gene counts have been stored in $final_output"
echo "$((($SECONDS-$checkpoint) / 60)) minutes and $((($SECONDS-$checkpoint) % 60)) seconds elapsed."
echo "Total time elapsed: $(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"

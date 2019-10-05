#!bin/bash

#this script will get data from RNA sequencing of two life stages of Tripanosoma bruscei bruscei,
#with three samples each, and count the mean count of reads aligning with each reference gene

#the 'fqfiles' file will guide a loop in order to process the sequencing pairs one at a time

guide="/localdisk/data/BPSM/Assignment1/fastq/fqfiles"
echo "$guide"


while read smpl_nmbr smpl_type pair1 pair2; do
 if test -f /localdisk/data/BPSM/Assignment1/fastq/$pair1; then
  echo "$pair1 exists"
  fi
 if test -f /localdisk/data/BPSM/Assignment1/fastq/$pair2; then
 echo "$pair2 exists"
  fi
done < $guide

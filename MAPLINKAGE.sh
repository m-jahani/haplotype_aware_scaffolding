#!/bin/bash

REF=$1      #/DATA/home/mjahani/TEST/AGA10_hap1.fasta
prime_3=$2  #/DATA/home/mjahani/test/3_primeBoundary.fasta
prime_5=$3  #/DATA/home/mjahani/test/5_primeBoundary.fasta
SAVE_DIR=$4 #/DATA/home/mjahani/test

cd $SAVE_DIR

bwa index -p $(basename ${REF%%.fasta}) -a bwtsw $REF
bwa mem -t 152 $(basename ${REF%%.fasta}) $prime_3 | samtools view -F 2308 | grep -v "^@" | awk '{print $1,"\t",$3,"\t",$4}' >${SAVE_DIR}/$(basename "${prime_3%%.fasta}")_$(basename "${REF%%.fasta}.sam")
bwa mem -t 152 $(basename ${REF%%.fasta}) $prime_5 | samtools view -F 2308 | grep -v "^@" | awk '{print $1,"\t",$3,"\t",$4}' >${SAVE_DIR}/$(basename "${prime_5%%.fasta}")_$(basename "${REF%%.fasta}.sam")
rm $(basename "${REF%%.fasta}")*

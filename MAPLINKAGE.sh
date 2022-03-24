#!/bin/bash

REF=/DATA/home/mjahani/TEST/AGA10_hap1.fasta
prime_3=/DATA/home/mjahani/test/3_primeBoundary.fasta
prime_5=/DATA/home/mjahani/test/5_primeBoundary.fasta
SAVE_DIR=/DATA/home/mjahani/test

cd $SAVE_DIR

bwa index -p $(basename ${REF%%.fasta}) -a bwtsw $REF
bwa mem -t 152 $(basename ${REF%%.fasta}) $prime_3 | samtools view -F 2308 | grep -v "^@" | awk '{print $1,"\t",$3,"\t",$4}' >${SAVE_DIR}/$(basename "${prime_3%%.fasta}.sam")
bwa mem -t 152 $(basename ${REF%%.fasta}) $prime_5 | samtools view -F 2308 | grep -v "^@" | awk '{print $1,"\t",$3,"\t",$4}' >${SAVE_DIR}/$(basename "${prime_5%%.fasta}.sam")

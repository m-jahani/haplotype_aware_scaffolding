#!/bin/bash

conda activate ASSEM_SCAFF2

ASSEMBLY=$1   #assembly name #Kalea.hic.hap1.p_ctg.fasta
PREFIX=$2     #scaffoldinf name in working directory #Kalea_hap1 in /home/mjahani/scratch/opt/juicer/scaffolding
JUCIER_DIR=$3 #/lustre04/scratch/mjahani/opt/juicer
HIC_R1=$4     #/DATA/home/mjahani/curation_AGA10/NS.1630.004.IDT_i7_195---IDT_i5_194.AGA10_HiC_R1.fastq.gz
HIC_R2=$5

[ -d ${JUCIER_DIR}/scaffolding/${PREFIX}/fastq ] || mkdir -p ${JUCIER_DIR}/scaffolding/${PREFIX}/fastq
cp -R -u -p $HIC_R1 ${JUCIER_DIR}/scaffolding/${PREFIX}/fastq
cp -R -u -p $HIC_R2 ${JUCIER_DIR}/scaffolding/${PREFIX}/fastq

cd ${JUCIER_DIR}/scaffolding/${PREFIX}

${JUCIER_DIR}/scripts/juicer.sh \
    -g $PREFIX \
    -s MboI \
    -t 152 \
    -z ${JUCIER_DIR}/references/${ASSEMBLY##*/} \
    -y ${JUCIER_DIR}/restriction_sites/${PREFIX}_MboI.txt \
    -p ${JUCIER_DIR}/chromosome_size/${PREFIX}.chrom.sizes \
    -D $JUCIER_DIR

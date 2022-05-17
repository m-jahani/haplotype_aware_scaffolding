#!/bin/bash

conda activate ASSEM_SCAFF2

ASSEMBLY=$1   #assembly name #Kalea.hic.hap1.p_ctg.fasta
PREFIX=$2     #scaffoldinf name in working directory #Kalea_hap1 in /home/mjahani/scratch/opt/juicer/scaffolding
JUCIER_DIR=$3 #/lustre04/scratch/mjahani/opt/juicer

[ -d ${JUCIER_DIR}/scaffolding/${PREFIX}/3D_DNA ] || mkdir -p ${JUCIER_DIR}/scaffolding/${PREFIX}/3D_DNA
cd ${JUCIER_DIR}/scaffolding/${PREFIX}/3D_DNA

bash ${JUCIER_DIR}/3d-dna/run-asm-pipeline.sh -i 5000 -r 0 --editor-repeat-coverage 5 --editor-coarse-resolution 100000 --editor-coarse-region 500000 ${JUCIER_DIR}/references/${ASSEMBLY##*/} ${JUCIER_DIR}/scaffolding/${PREFIX}/aligned/merged_nodups.txt

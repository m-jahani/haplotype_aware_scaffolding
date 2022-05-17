#!/bin/bash

conda activate ASSEM_SCAFF2

ASSEMBLY=$1   #assembly name #Kalea.hic.hap1.p_ctg.fasta
PREFIX=$2     #scaffoldinf name in working directory #Kalea_hap1 in /home/mjahani/scratch/opt/juicer/scaffolding
JUCIER_DIR=$3 #/lustre04/scratch/mjahani/opt/juicer

cp -R -u -p $ASSEMBLY ${JUCIER_DIR}/references/
bwa index ${JUCIER_DIR}/references/${ASSEMBLY##*/}

cd ${JUCIER_DIR}/restriction_sites
python ${JUCIER_DIR}/misc/generate_site_positions.py MboI $PREFIX ${JUCIER_DIR}/references/${ASSEMBLY##*/}
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${PREFIX}_MboI.txt >${JUCIER_DIR}/chromosome_size/${PREFIX}.chrom.sizes

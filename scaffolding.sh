#!/bin/bash
bin=$1 #/DATA/home/mjahani/CURATION/bin #directory of the scripts #/DATA/home/mjahani/CURATION/bin
ID=$2 #AGA10_hap1 #Id for the sample #AGA10_hap1
ASSEMBLY=$3 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap1.p_ctg.fasta #AGA10.hic.hap1.p_ctg.fasta
JUCIER_DIR=$4 #/DATA/home/mjahani/opt/juicer #/DATA/home/mjahani/opt/juicer
HIC_R1=$5 #/DATA/home/mjahani/ASSEMBLIES/AGA10/HiC/NS.1630.004.IDT_i7_195---IDT_i5_194.AGA10_HiC_R1.fastq.gz
HIC_R2=$6 #/DATA/home/mjahani/ASSEMBLIES/AGA10/HiC/NS.1630.004.IDT_i7_195---IDT_i5_194.AGA10_HiC_R1.fastq.gz
MAPQ=$7 # for HiC alignments in 3D_DNA
RESULT_DIR=$8
#preperation for first round scaffolding: indexes the assembly, finds the restriction enzyme locations, calculate the contig sizes
bash ${bin}/JUICER_PREP.sh $ASSEMBLY $ID $JUCIER_DIR ${RESULT_DIR}

#Script to map HiC data on assembly for 3D DNA scaffolding
bash ${bin}/JUICER.sh $ASSEMBLY $ID $JUCIER_DIR $HIC_R1 $HIC_R2 ${RESULT_DIR}
#First round of 3D DNA scaffolding
bash ${bin}/3DDNA.sh $ASSEMBLY $ID $JUCIER_DIR ${MAPQ} ${RESULT_DIR}

#!/bin/bash
bin=/DATA/home/mjahani/CURATION/bin #directory of the scripts #/DATA/home/mjahani/CURATION/bin
SAMPLE=$1                           #3AGA10
ASSEMBLY_hap1=$2                    #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap1.p_ctg.fasta #AGA10.hic.hap1.p_ctg.fasta
ASSEMBLY_hap2=$3                    #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap2.p_ctg.fasta
RESULT_DIR=$4
CURATED_ASSEMBLY=$5

#mix assembly to fasta
[ -d ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP

Rscript ${bin}/ASM2FASTA_GAP.R \
    $CURATED_ASSEMBLY \
    ${SAMPLE}_FINAL_GAP \
    ${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP

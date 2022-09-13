#!/bin/bash
bin=/DATA/home/mjahani/CURATION/bin #directory of the scripts #/DATA/home/mjahani/CURATION/bin
SAMPLE=$1                           #3AGA10
ASSEMBLY_hap1=$2                    #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap1.p_ctg.fasta #AGA10.hic.hap1.p_ctg.fasta
ASSEMBLY_hap2=$3                    #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap2.p_ctg.fasta
RESULT_DIR=$4
CURATED_ASSEMBLY=$5

CS10_ASSEMBLY=${RESULT_DIR}/CS10_assembly/GCF_900626175.2_cs10_10CHR.fasta

#mix assembly to fasta
[ -d ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP

Rscript ${bin}/ASM2FASTA_GAP.R \
    $CURATED_ASSEMBLY \
    ${SAMPLE}_FINAL_GAP \
    ${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP

GENOME_ASSEMBLY_HAP1=${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP/${SAMPLE}_FINAL_GAP_hap1.reviewed.chr_assembled.fasta
GENOME_ASSEMBLY_HAP2=${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_GAP/${SAMPLE}_FINAL_GAP_hap2.reviewed.chr_assembled.fasta
#####Hap1
bash ${bin}/minimap.sh \
    $GENOME_ASSEMBLY_HAP1 \
    $CS10_ASSEMBLY \
    ${RESULT_DIR}/${SAMPLE}/minimap

[ -d ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_ASSEMBLY ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_ASSEMBLY
Rscript ${bin}/CHRID.R \
    ${RESULT_DIR}/${SAMPLE}/minimap/$(basename "${GENOME_ASSEMBLY_HAP1%%.fasta}")_$(basename "${CS10_ASSEMBLY%%.fasta}").paf \
    $GENOME_ASSEMBLY_HAP1 \
    ${SAMPLE}_Hap1 \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_ASSEMBLY

#####Hap2
bash ${bin}/minimap.sh \
    $GENOME_ASSEMBLY_HAP2 \
    $CS10_ASSEMBLY \
    ${RESULT_DIR}/${SAMPLE}/minimap

Rscript ${bin}/CHRID.R \
    ${RESULT_DIR}/${SAMPLE}/minimap/$(basename "${GENOME_ASSEMBLY_HAP2%%.fasta}")_$(basename "${CS10_ASSEMBLY%%.fasta}").paf \
    $GENOME_ASSEMBLY_HAP2 \
    ${SAMPLE}_Hap2 \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/FINAL_ASSEMBLY

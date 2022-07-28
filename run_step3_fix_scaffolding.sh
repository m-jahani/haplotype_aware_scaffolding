#!/bin/bash

bin=/DATA/home/mjahani/CURATION/bin #directory of the scripts #/DATA/home/mjahani/CURATION/bin
ID_HAP1=$1                          #AGA10_hap1 #ID for the sample #AGA10_hap1
ID_HAP2=$2                          #AGA10_hap2 #ID for the sample #AGA10_hap2
SAMPLE=$3                           #3AGA10
ASSEMBLY_hap1=$4                    #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap1.p_ctg.fasta #AGA10.hic.hap1.p_ctg.fasta
ASSEMBLY_hap2=$5                    #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap2.p_ctg.fasta
RESULT_DIR=$6                       #/DATA/home/mjahani/CURATION
curation_round=$7
CURATED_ASSEMBLY=$8 #/DATA/home/mjahani/CURATION/AGA10/HiC_map/AGA10.hic.hap1.p_ctg_AGA10.hic.hap2.p_ctg.1.review.assembly

TELOMERE_WINSIZE=200000
#mix assembly to fasta
[ -d ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round} ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}
Rscript ${bin}/ASM2FASTA.R \
    $CURATED_ASSEMBLY \
    ${SAMPLE}_Round${curation_round} \
    ${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}

#SYNTENY
bash ${bin}/minimap.sh \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}/${SAMPLE}_Round${curation_round}_hap1.reviewed.chr_assembled.fasta \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}/${SAMPLE}_Round${curation_round}_hap2.reviewed.chr_assembled.fasta \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}

#drawing plots
[ -d ${RESULT_DIR}/${SAMPLE}/PLOT/Round${curation_round} ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/PLOT/Round${curation_round}

Rscript ${bin}/interactive_plot.R \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original/${SAMPLE}_hap12.reviwed_contig_chr_coord \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}/${SAMPLE}_Round${curation_round}_hap12.reviwed_contig_chr_coord \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled.GeneticMap \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}.reviewed.chr_assembled.GeneticMap \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled.recombination \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}.reviewed.chr_assembled.recombination \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}_telomeric_repeat_${TELOMERE_WINSIZE}windows.bed \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}_telomeric_repeat_${TELOMERE_WINSIZE}windows.bed \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled.EDTA.bed \
    ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}.reviewed.chr_assembled.EDTA.bed \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}/${SAMPLE}_Round${curation_round}_hap1.reviewed.chr_assembled_${SAMPLE}_Round${curation_round}_hap2.reviewed.chr_assembled.paf \
    ${RESULT_DIR}/${SAMPLE}/PLOT/Round${curation_round}

Rscript ${bin}/breakingPoint.R \
    ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/Round${curation_round}/${SAMPLE}_Round${curation_round}_hap1.reviewed.chr_assembled_${SAMPLE}_Round${curation_round}_hap2.reviewed.chr_assembled.paf \
    ${RESULT_DIR}/${SAMPLE}/PLOT/Round${curation_round}

#!/bin/bash
bin=/DATA/home/mjahani/CURATION/bin #directory of the scripts #/DATA/home/mjahani/CURATION/bin
ID_HAP1=$1 #AGA10_hap1 #ID for the sample #AGA10_hap1
ID_HAP2=$2 #AGA10_hap2 #ID for the sample #AGA10_hap2
SAMPLE=$3 #AGA10
ASSEMBLY_hap1=$4 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap1.p_ctg.fasta #AGA10.hic.hap1.p_ctg.fasta
ASSEMBLY_hap2=$5 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap2.p_ctg.fasta
RESULT_DIR=$6 #/DATA/home/mjahani/CURATION
PRIME_3=$7 #/DATA/home/mjahani/CURATION/GeneticMap/3_primeBoundary.fasta
PRIME_5=$8 #/DATA/home/mjahani/CURATION/GeneticMap/5_primeBoundary.fasta
LINKAGE_DATA=$9 #/DATA/home/mjahani/CURATION/GeneticMap/MSTmap_4_GAPP_project.csv

TELOMERE_WINSIZE=200000
REVIEWED_FASTA_HAP1=${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original/${ID_HAP1}.reviewed.chr_assembled.fasta
REVIEWED_FASTA_HAP2=${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original/${ID_HAP2}.reviewed.chr_assembled.fasta

#mix assembly files for JUCIE BOX
Rscript ${bin}/mix_reviewed_assemblies.R \
${RESULT_DIR}/${SAMPLE}/HiC_map/$(basename ${ASSEMBLY_hap1%%.fasta}).0.review.assembly \
${RESULT_DIR}/${SAMPLE}/HiC_map/$(basename ${ASSEMBLY_hap2%%.fasta}).0.review.assembly \
${RESULT_DIR}/${SAMPLE}/minimap/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).paf \
${RESULT_DIR}/${SAMPLE}/HiC_map 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_${ID_HAP2}_mix_reviewed_assemblies_original.LOG

#mix assembly to fasta
[ -d ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original
Rscript ${bin}/ASM2FASTA.R ${RESULT_DIR}/${SAMPLE}/HiC_map/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).review.assembly $SAMPLE ${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta ${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_${ID_HAP2}_ASM2FASTA_original.LOG

[ -d ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO
Linkage map and recombination hap1
bash ${bin}/MAPLINKAGE.sh $REVIEWED_FASTA_HAP1 $PRIME_3 $PRIME_5 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_geneticmap_maping.LOG

Rscript ${bin}/recombination.R \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_3%%.fasta}")_${ID_HAP1}.reviewed.chr_assembled.sam \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_5%%.fasta}")_${ID_HAP1}.reviewed.chr_assembled.sam \
$LINKAGE_DATA \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_geneticmap_recombination.LOG

rm ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_3%%.fasta}")_${ID_HAP1}.reviewed.chr_assembled.sam
rm ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_5%%.fasta}")_${ID_HAP1}.reviewed.chr_assembled.sam

#Linkage map and recombination hap2
bash ${bin}/MAPLINKAGE.sh $REVIEWED_FASTA_HAP2 $PRIME_3 $PRIME_5 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP2}_geneticmap_maping.LOG

Rscript ${bin}/recombination.R \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_3%%.fasta}")_${ID_HAP2}.reviewed.chr_assembled.sam \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_5%%.fasta}")_${ID_HAP2}.reviewed.chr_assembled.sam \
$LINKAGE_DATA \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO  2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP2}_geneticmap_recombination.LOG

rm ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_3%%.fasta}")_${ID_HAP2}.reviewed.chr_assembled.sam
rm ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/$(basename "${PRIME_5%%.fasta}")_${ID_HAP2}.reviewed.chr_assembled.sam

#Telomere hap1
bash ${bin}/TELOMERE.sh $REVIEWED_FASTA_HAP1 $TELOMERE_WINSIZE $ID_HAP1 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO  2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_telomere.LOG
#Telomere hap1
bash ${bin}/TELOMERE.sh $REVIEWED_FASTA_HAP2 $TELOMERE_WINSIZE $ID_HAP2 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP2}_telomere.LOG

#EDTA hap1
bash ${bin}/EDTA.sh $REVIEWED_FASTA_HAP1 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_EDTA.LOG &
#EDTA hap2
bash ${bin}/EDTA.sh $REVIEWED_FASTA_HAP2 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO 2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP2}_EDTA.LOG

#SYNTENY
bash ${bin}/minimap.sh $REVIEWED_FASTA_HAP1 $REVIEWED_FASTA_HAP2 ${RESULT_DIR}/${SAMPLE}/CONTIG_INFO  2> ${RESULT_DIR}/${SAMPLE}/LOG/${ID_HAP1}_${ID_HAP2}_stnteny_minimap.LOG

[ -d ${RESULT_DIR}/${SAMPLE}/PLOT/original ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/PLOT/original
#drawing plots
Rscript ${bin}/interactive_plot.R \
${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original/${SAMPLE}_hap12.reviwed_contig_chr_coord \
${RESULT_DIR}/${SAMPLE}/ASM2FASTA/original/${SAMPLE}_hap12.reviwed_contig_chr_coord \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled.GeneticMap \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}.reviewed.chr_assembled.GeneticMap \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled.recombination \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}.reviewed.chr_assembled.recombination \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}_telomeric_repeat_${TELOMERE_WINSIZE}windows.bed \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}_telomeric_repeat_${TELOMERE_WINSIZE}windows.bed \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled.EDTA.bed \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP2}.reviewed.chr_assembled.EDTA.bed \
${RESULT_DIR}/${SAMPLE}/CONTIG_INFO/${ID_HAP1}.reviewed.chr_assembled_${ID_HAP2}.reviewed.chr_assembled.paf \
${RESULT_DIR}/${SAMPLE}/PLOT/original

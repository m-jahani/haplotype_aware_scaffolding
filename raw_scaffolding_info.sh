#!/bin/bash

#raw assembly
RAW_SCAFFOLDS_HAP1=$1 #/DATA/home/mjahani/curation_AGA10/AGA10hap1.reviewed.chr_assembled.fasta
RAW_SCAFFOLDS_HAP2=$2 #/DATA/home/mjahani/curation_AGA10/AGA10hap1.reviewed.chr_assembled.fasta
#assembly graph
RAW_CONTIGS_HAP1_GRAPH=$3 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_HIFIasm_ASSEMBLY/AGA10.hic.hap1.p_ctg.gfa
RAW_CONTIGS_HAP2_GRAPH=$4 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_HIFIasm_ASSEMBLY/AGA10.hic.hap2.p_ctg.gfa
#linkage map markers
prime_3=$5      #/DATA/home/mjahani/test/3_primeBoundary.fasta
prime_5=$6      #/DATA/home/mjahani/test/5_primeBoundary.fasta
LINKAGE_DATA=$7 #/DATA/home/mjahani/curation_AGA10/MSTmap_4_GAPP_project.csv
SAVE_DIR=$8     #/DATA/home/mjahani/curation_AGA10/TEMP_TEST/result

[ -d ${SAVE_DIR}/LINKAGE_MAP ] || mkdir -p ${SAVE_DIR}/LINKAGE_MAP
bash MAPLINKAGE.sh $RAW_SCAFFOLDS_HAP1 $prime_3 $prime_5 ${SAVE_DIR}/LINKAGE_MAP
Rscript recombination.R ${SAVE_DIR}/LINKAGE_MAP/$(basename "${prime_3%%.fasta}")_$(basename "${RAW_SCAFFOLDS_HAP1%%.fasta}.sam") ${SAVE_DIR}/LINKAGE_MAP/$(basename "${prime_5%%.fasta}")_$(basename "${RAW_SCAFFOLDS_HAP1%%.fasta}.sam") $LINKAGE_DATA ${SAVE_DIR}/LINKAGE_MAP

bash MAPLINKAGE.sh $RAW_SCAFFOLDS_HAP2 $prime_3 $prime_5 ${SAVE_DIR}/LINKAGE_MAP
Rscript recombination.R ${SAVE_DIR}/LINKAGE_MAP/$(basename "${prime_3%%.fasta}")_$(basename "${RAW_SCAFFOLDS_HAP2%%.fasta}.sam") ${SAVE_DIR}/LINKAGE_MAP/$(basename "${prime_5%%.fasta}")_$(basename "${RAW_SCAFFOLDS_HAP2%%.fasta}.sam") $LINKAGE_DATA ${SAVE_DIR}/LINKAGE_MAP

[ -d ${SAVE_DIR}/TELOMERE ] || mkdir -p ${SAVE_DIR}/TELOMERE
bash TELOMERE.sh $RAW_SCAFFOLDS_HAP1 200000 $(basename ${RAW_SCAFFOLDS_HAP1%.fasta}) ${SAVE_DIR}/TELOMERE
bash TELOMERE.sh $RAW_SCAFFOLDS_HAP2 200000 $(basename ${RAW_SCAFFOLDS_HAP2%.fasta}) ${SAVE_DIR}/TELOMERE

[ -d ${SAVE_DIR}/DEPTH ] || mkdir -p ${SAVE_DIR}/DEPTH
bash GFA2DEPTH.sh $RAW_CONTIGS_HAP1_GRAPH ${SAVE_DIR}/DEPTH
bash GFA2DEPTH.sh $RAW_CONTIGS_HAP2_GRAPH ${SAVE_DIR}/DEPTH

[ -d ${SAVE_DIR}/EDTA_LTR ] || mkdir -p ${SAVE_DIR}/EDTA_LTR
bash EDTA.sh $RAW_SCAFFOLDS_HAP1 ${SAVE_DIR}/EDTA_LTR &
bash EDTA.sh $RAW_SCAFFOLDS_HAP2 ${SAVE_DIR}/EDTA_LTR

[ -d ${SAVE_DIR}/SYNTENY ] || mkdir -p ${SAVE_DIR}/SYNTENY
bash minimap.sh $RAW_SCAFFOLDS_HAP1 $RAW_SCAFFOLDS_HAP2 ${SAVE_DIR}/SYNTENY

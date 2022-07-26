#!/bin/bash
bin=/DATA/home/mjahani/CURATION/bin #directory of the scripts #/DATA/home/mjahani/CURATION/bin
ID_HAP1=$1 #AGA10_hap1 #ID for the sample #AGA10_hap1
ID_HAP2=$2 #AGA10_hap2 #ID for the sample #AGA10_hap2
SAMPLE=$3 #AGA10
ASSEMBLY_hap1=$4 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap1.p_ctg.fasta #AGA10.hic.hap1.p_ctg.fasta
ASSEMBLY_hap2=$5 #/DATA/home/mjahani/ASSEMBLIES/AGA10/AGA10_ASSEMBLY_RESULT_STAT/AGA10.hic.hap2.p_ctg.fasta
JUCIER_DIR=$6 #/DATA/home/mjahani/opt/juicer #/DATA/home/mjahani/opt/juicer
RESULT_DIR=$7 #/DATA/home/mjahani/CURATION
HIC_R1=$8 #/DATA/home/mjahani/ASSEMBLIES/AGA10/HiC/NS.1630.004.IDT_i7_195---IDT_i5_194.AGA10_HiC_R1.fastq.gz
HIC_R2=$9 #/DATA/home/mjahani/ASSEMBLIES/AGA10/HiC/NS.1630.004.IDT_i7_195---IDT_i5_194.AGA10_HiC_R2.fastq.gz


[ -d ${RESULT_DIR}/${SAMPLE} ] || mkdir -p ${RESULT_DIR}/${SAMPLE}
RESULT_DIR_ID=${RESULT_DIR}/${SAMPLE}
#hap1
bash ${bin}/scaffolding.sh $bin $ID_HAP1 $ASSEMBLY_hap1 $JUCIER_DIR $HIC_R1 $HIC_R2 1 $RESULT_DIR_ID &
#hap2
bash ${bin}/scaffolding.sh $bin $ID_HAP2 $ASSEMBLY_hap2 $JUCIER_DIR $HIC_R1 $HIC_R2 1 $RESULT_DIR_ID &
#mix hap
[ -d ${RESULT_DIR}/${SAMPLE}/sequences ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/sequences
bash ${bin}/cat_fasta.sh $ASSEMBLY_hap1 $ASSEMBLY_hap2 ${RESULT_DIR}/${SAMPLE}/sequences
bash ${bin}/scaffolding.sh $bin ${ID_HAP1}_${ID_HAP2} ${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta $JUCIER_DIR $HIC_R1 $HIC_R2 10 $RESULT_DIR_ID
#mini map alignment
[ -d ${RESULT_DIR}/${SAMPLE}/minimap ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/minimap
bash ${bin}/minimap.sh $ASSEMBLY_hap1 $ASSEMBLY_hap2 ${RESULT_DIR}/${SAMPLE}/minimap

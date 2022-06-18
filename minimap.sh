#!/bin/bash

conda activate scaff_env

ASSEMBLY1=$1 #/DATA/home/mjahani/opt/juicer/scaffolding/AGA10_hap1/3D_DNA/AGA10.hic.hap1.p_ctg.final.fasta
ASSEMBLY2=$2 #/DATA/home/mjahani/opt/juicer/scaffolding/AGA10_hap2/3D_DNA/AGA10.hic.hap2.p_ctg.final.fasta
SAVE_DIR=$3  #/DATA/home/mjahani/curation_AGA10
ASSEMBLY1_FILE=${ASSEMBLY1##*/}
ASSEMBLY2_FILE=${ASSEMBLY2##*/}

/DATA/home/mjahani/bin/minimap2/minimap2 -x asm10 -t 152 $ASSEMBLY1 $ASSEMBLY2 >${SAVE_DIR}/$(basename "${RAW_CONTIGS_HAP1%%.fasta}")_$(basename "${RAW_CONTIGS_HAP2%%.fasta}").paf

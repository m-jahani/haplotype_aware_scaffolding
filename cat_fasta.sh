#!/bin/bash
FASTA1=$1
FASTA2=$2
SAVE_DIR=$3
FASTA1_FILE=${FASTA1##*/}
FASTA2_FILE=${FASTA2##*/}

cat $FASTA1 $FASTA2 >${SAVE_DIR}/${FASTA1_FILE%%.fasta}_${FASTA2_FILE%%.fasta}.fasta

#!/bin/bash

source activate scaff

FASTA=$1
WIN_SIZE=$2 #tested on 20000
PREFIX=$3
SAVE_DIR=$4

tidk search --fasta $FASTA --string TTTAGGG --window $WIN_SIZE --output $PREFIX --dir $SAVE_DIR
awk -F',' -v WINSIZE=${WIN_SIZE} '{print $1, $2-WINSIZE, $2, $4 }' ${SAVE_DIR}/${PREFIX}_telomeric_repeat_windows.csv | sed 1d >${SAVE_DIR}/${PREFIX}_telomeric_repeat_${WIN_SIZE}windows.bed

#!/bin/bash
conda activate ASSEM_SCAFF2
FASTA=$1
SAVEDIR=$2

FASTA_BASE=${FASTA##*/}
grep "^>" $FASTA | sed 's/^>//g' >${SAVEDIR}/${FASTA_BASE%%.fasta}_CHR_LISR

while read SCAFFOLD; do
    echo $SCAFFOLD >${SAVEDIR}/${SCAFFOLD}.txt
    seqtk subseq $FASTA ${SAVEDIR}/${SCAFFOLD}.txt >${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa
    rm ${SAVEDIR}/${SCAFFOLD}.txt
    trf ${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa 2 5 7 80 10 50 2000 -h
    awk -v i="$SCAFFOLD" '{print i,"\t",$1,"\t",$2}' ${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa.2.5.7.80.10.50.2000.dat | awk '{ if ($2 == ($2+0)) print $0 }' | sed 's/ //g' >${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.bed
    rm ${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa ${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa.2.5.7.80.10.50.2000.dat
    cat ${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.bed >>${SAVEDIR}/${FASTA_BASE%%.fasta}.bed
    rm ${SAVEDIR}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.bed
done <${SAVEDIR}/${FASTA_BASE%%.fasta}_CHR_LISR

bedtools maskfasta -soft -fi $FASTA -bed ${SAVEDIR}/${FASTA_BASE%%.fasta}.bed -fo ${SAVEDIR}/${FASTA_BASE%%.fasta}.softmask.fatsa
rm ${SAVEDIR}/${FASTA_BASE%%.fasta}.bed ${SAVEDIR}/${FASTA_BASE%%.fasta}_CHR_LISR

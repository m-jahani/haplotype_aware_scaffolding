#!/bin/bash
conda activate ASSEM_SCAFF2
FASTA=$1

grep "^>" $FASTA | sed 's/^>//g' >${FASTA%%.fasta}_CHR_LISR

while read SCAFFOLD; do
    echo $SCAFFOLD >${SCAFFOLD}.txt
    seqtk subseq $FASTA ${SCAFFOLD}.txt >${FASTA%%.fasta}_${SCAFFOLD}.fatsa
    rm ${SCAFFOLD}.txt
    trf ${FASTA%%.fasta}_${SCAFFOLD}.fatsa 2 5 7 80 10 50 2000 -h
    awk -v i="$SCAFFOLD" '{print i,"\t",$1,"\t",$2}' ${FASTA%%.fasta}_${SCAFFOLD}.fatsa.2.5.7.80.10.50.2000.dat | awk '{ if ($2 == ($2+0)) print $0 }' | sed 's/ //g' >${FASTA%%.fasta}_${SCAFFOLD}.bed
    rm ${FASTA%%.fasta}_${SCAFFOLD}.fatsa ${FASTA%%.fasta}_${SCAFFOLD}.fatsa.2.5.7.80.10.50.2000.dat
    cat ${FASTA%%.fasta}_${SCAFFOLD}.bed >>${FASTA%%.fasta}.bed
    rm ${FASTA%%.fasta}_${SCAFFOLD}.bed
done <${FASTA%%.fasta}_CHR_LISR

bedtools maskfasta -soft -fi $FASTA -bed ${FASTA%%.fasta}.bed -fo ${FASTA%%.fasta}.softmask.fatsa
rm ${FASTA%%.fasta}.bed ${FASTA%%.fasta}_CHR_LISR

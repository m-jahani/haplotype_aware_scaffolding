#!/bin/bash
eval "$(command conda 'shell.bash' 'hook' 2>/dev/null)"
conda activate EDTA
FASTA=$1 #/DATA/home/mjahani/curation_AGA10/AGA10hap1.reviewed.chr_assembled.fasta
SAVEDIR=$2

FASTA_BASE=${FASTA##*/}

bioawk -c fastx '{ print $name, length($seq) }' $FASTA | sort -k2nr | cut -f1 | head -10 >${SAVEDIR}/${FASTA_BASE%%.fasta}_CHR_LIST

while read SCAFFOLD; do
    echo $SCAFFOLD >${SAVEDIR}/${SCAFFOLD}.txt
    mkdir -p ${SAVEDIR}/${FASTA_BASE%%.fasta}/${SCAFFOLD}
    seqtk subseq $FASTA ${SAVEDIR}/${SCAFFOLD}.txt >${SAVEDIR}/${FASTA_BASE%%.fasta}/${SCAFFOLD}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa
    rm ${SAVEDIR}/${SCAFFOLD}.txt
    cd ${SAVEDIR}/${FASTA_BASE%%.fasta}/${SCAFFOLD}
    singularity exec /DATA/home/mjahani/bin/EDTA/EDTA.sif EDTA.pl --overwrite 0 --genome ${SAVEDIR}/${FASTA_BASE%%.fasta}/${SCAFFOLD}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa --sensitive 0 --anno 1 --evaluate 0 --threads 76
    awk -v var=$SCAFFOLD '{print var,"\t",$3,"\t",$4,"\t",$5}' ${SAVEDIR}/${FASTA_BASE%%.fasta}/${SCAFFOLD}/${FASTA_BASE%%.fasta}_${SCAFFOLD}.fatsa.mod.EDTA.TEanno.gff3 >>${SAVEDIR}/${FASTA_BASE%%.fasta}_10chr.bed
done <${SAVEDIR}/${FASTA_BASE%%.fasta}_CHR_LIST

rm ${SAVEDIR}/${FASTA_BASE%%.fasta}_CHR_LIST

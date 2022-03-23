#!/bin/bash
GFA=$1

awk '/^A/{print $2,"\t",$3,"\t",$3+$7,"\t",$8}' $GFA | sed 's/ //g' >${GFA}_READ.bed
bedtools sort -i ${GFA}_READ.bed >${GFA}_READ.sorted.bed
rm ${GFA}_READ.bed
awk '/^S/{print $2,"\t",$4}' $GFA | sed 's/LN:i://g' | sed 's/ //g' >${GFA}_contigs
bedtools genomecov -i ${GFA}_READ.sorted.bed -g ${GFA}_contigs -bga >${GFA%%.gfa}.depth
rm ${GFA}_READ.sorted.bed ${GFA}_contigs

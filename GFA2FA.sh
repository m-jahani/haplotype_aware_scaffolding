#!/bin/bash

GFA=$1
awk '/^S/{print ">"$2;print $3}' $GFA >${GFA%%gfa}fasta

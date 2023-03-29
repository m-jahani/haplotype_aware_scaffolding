# haplotype_aware_scaffolding

a series of custom scripts for haplotype aware manual curation of genome scaffolds

## Description

This pipeline was designed for haplotype-aware scaffolding of Cannabis genome assemblies using Juicer and 3DDNA for mapping Hi-C and scaffolding the genome, respectively. Unlike traditional curation of one haplotype at a time, this pipeline generates inputs for Juicebox to view two haplotypes simultaneously. Additionally, it generates interactive plots for both haplotypes to visualize genetic distance of markers, recombination rate, telomere signal frequency, LTR frequency, and alignment plot of the two haplotypes. These plots, alongside Hi-C plots, provide enough information for manual curation of Cannabis genome assemblies in Juicebox. After one or several rounds of curation, the pipeline can finalize genomes by adding NNNN space between contigs and assigning chromosome numbers to scaffolds based on a pre-defined genome assembly (in this case, CS10).

#                                           Under Development

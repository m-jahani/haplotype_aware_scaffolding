# Haplotype Aware Scaffolding

A series of custom scripts for haplotype aware manual curation of Cannabis genome assemblies

## Description

This pipeline was designed for haplotype-aware scaffolding of Cannabis genome assemblies using Juicer and 3DDNA for mapping Hi-C and scaffolding the genome, respectively. Unlike traditional curation of one haplotype at a time, this pipeline generates inputs for Juicebox to view two haplotypes simultaneously. Additionally, it generates interactive plots for both haplotypes to visualize genetic distance of markers, recombination rate, telomere signal frequency, LTR frequency, and alignment plot of the two haplotypes. These plots, alongside Hi-C plots, provide enough information for manual curation of Cannabis genome assemblies in Juicebox. After one or several rounds of curation, the pipeline can finalize genomes by adding NNNN space between contigs and assigning chromosome numbers to scaffolds based on a pre-defined genome assembly (in this case, CS10).


## Step 1. HiC mapping and Scaffolding

This step focuses on mapping Hi-C reads to Haplotype 1 and Haplotype 2 of the Cannabis genome assembly separately, and then runs contig scaffolding. Additionally, it combines both haplotypes in a single FASTA file and repeats the mapping and scaffolding step. To identify corresponding scaffolds in the haplotypes of a genome assembly, the pipeline aligns the haplotypes against each other using minimap2.

It is important to note that, in 3D-DNA scaffolding, the parameter `-r 0` was used. `r` defines the number of misjoin correction rounds, and by setting it to `0`, we skipped misjoin correction since 3D-DNA tends to over-correct plant genome assemblies.

Step 1 can be executed with [run_step1_scaffolding.sh](https://github.com/m-jahani/haplotype_aware_scaffolding/blob/main/run_step1_scaffolding.sh) as follow:
```
bash run_step1_scaffolding.sh \
ID_HAP1 \ 
ID_HAP2 \ 
SAMPLE \
ASSEMBLY_hap1 \
ASSEMBLY_hap2 \
JUCIER_DIR \
RESULT_DIR \
HIC_R1 \
HIC_R2
```

Where each variable can be defined as follow:
```
ID_HAP1 = ID for haplotype one, example: MK_ultra_hap1
ID_HAP2 = ID for haplotype two, example: MK_ultra_hap2
SAMPLE = ID for the genotype, example: MK_ultra
ASSEMBLY_hap1 = genome assembly of haplotype 1 in fasta format, example: MK_ultra.hic.hap1.p_ctg.fasta
ASSEMBLY_hap2 = genome assembly of haplotype 2 
in fasta format, example: MK_ultra.hic.hap2.p_ctg.fasta
JUCIER_DIR = Path to where Juicer was installed
RESULT_DIR = path for saving results
HIC_R1 = HiC read `R1` in fasta format
HIC_R2 = HiC read `R2` in fasta format
```

## Step 2. Drawing assembly curation guide plots
When the first step of genome scaffolding is done, we have to  draw curation guide plots (interactive information plots on haplotypes, and two haplotype HiC heatmaps) to get the enough information for manual curations. This step of pipline runs its different modules to capture this information and put them in plot format.

### Step 2.1
Generating a heatmap of both haplotypes in a single plot for Juicebox can be achieved using a custom R script (mix_reviewed_assemblies.R)[https://github.com/m-jahani/haplotype_aware_scaffolding/blob/main/mix_reviewed_assemblies.R]. This script takes the outputs of the first step from '3D-DNA' and 'minimap2' and combines the contig information data of the two haplotypes to generate a heatmap for 'Juicebox' with both haplotypes simultaneously in a plot.

### Step 2.2
Use the custom R script [ASM2FASTA.R](https://github.com/m-jahani/haplotype_aware_scaffolding/blob/main/ASM2FASTA.R) to convert the output of Step 2.1 to fasta format efficiently. The output of 3DDNA with the "*assembly" suffix is essentially a text file containing information on the order and orientation of contigs. Since this file can be edited in Juicebox, we can always obtain the fasta file after any curation.

### Step 2.3
[MAPLINKAGE.sh](https://github.com/m-jahani/haplotype_aware_scaffolding/blob/main/MAPLINKAGE.sh) calculates position of each genetic marker (from the genetic map) on each haplotype of the assembly. 

#                                           Under Development

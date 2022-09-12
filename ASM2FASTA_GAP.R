library(tidyverse)
library(tidysq)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

MIX_ASSEM <- args[1] #"/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/AGA10.hic.hap1_AGA10.hic.hap2.review.assembly"
PREFIX <- args[2] #AGA10
MIX_FASTA <- args[3]
SAVE_DIR <- args[4]

# MIX_ASSEM <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/MK-ULTRA/HiC/MK_ultra.hic.hap1.p_ctg_MK_ultra.hic.hap2.p_ctg.2.review.assembly"
# PREFIX <- "MK_ultra_Final"
# MIX_FASTA <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/MK-ULTRA/MK_ultra.hic.hap1.p_ctg_MK_ultra.hic.hap2.p_ctg.fasta"
#SAVE_DIR <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/MK-ULTRA"
####################################################################Read Data#####################################################################
# read *reviewed assembly hap1_hap2 ()assembly file after juicebox curation
read.table(MIX_ASSEM,
      fill = TRUE,
      col.names = paste("V", 1:100, sep = "")# with the assumption that each chromosome does not contain more than 100 conigs
) -> HAP12_assembly

read_fasta(MIX_FASTA) -> HAP12_FASTA #merged hap1 and hap2 contigs
##################################################Contigs/fragment/Debris order in assembly file####################################################
#fragmented contigs
HAP12_assembly %>%
  filter(grepl("^>",V1)) %>%
  select(fragment_ID = V1,
         original_order = V2,
         length = V3)   %>%
  separate(fragment_ID,
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>%
  filter(grepl("fragment",lable1)) ->  FRAGMENTS

#Run if some of the contigs are fragmented
if (nrow(FRAGMENTS) > 0 ) {

#calculate the coordination of the fragmented contigs in the original contig
FRAGMENTS %>%
  separate(lable1,
           into = c("Fragment","Fragment_order")) %>%
  group_by(contig) %>%
  arrange(as.numeric(Fragment_order)) %>%
  mutate(end=cumsum(length)) %>%
  ungroup() %>%
  #mutate(end = ifelse(Fragment_order == 1, length, lag(length) + length)) %>%
  mutate(start = (end-length) + 1) %>%
  mutate(fragment_ID = gsub(">","",fragment_ID)) %>%
  mutate(contig = gsub(">","",contig)) %>%
  select(contig,
         fragment_ID,
         original_order,
         start,
         end,
         length) -> coord_fragments

#extract hw sequence of the fragmented contigs from the original contig
FRAGMENTED_CONTIGS <- NULL
for (CONTIG in 1:nrow(coord_fragments)) {
  HAP12_FASTA %>%
    filter(name == as.character(coord_fragments[CONTIG,1])) %>%
    pull(sq) %>%
    bite(as.numeric(coord_fragments[CONTIG,4]):as.numeric(coord_fragments[CONTIG,5])) %>%
    enframe() %>%
    mutate(name = as.character(coord_fragments[CONTIG,2])) %>%
    select(sq = value,
           name) %>%
    bind_rows(.,
          FRAGMENTED_CONTIGS) -> FRAGMENTED_CONTIGS
}

#remove the the original contig and replace the fragments
HAP12_FASTA %>%
  filter(!name %in%  pull(distinct(coord_fragments,contig))) %>%
  bind_rows(.,
            FRAGMENTED_CONTIGS) -> HAP12_FASTA_FRAG

      rm(HAP12_FASTA,
        coord_fragments,
        FRAGMENTED_CONTIGS,
        CONTIG)
      HAP12_FASTA_FRAG -> HAP12_FASTA
      rm(HAP12_FASTA_FRAG)
}else {
  print("NO fragmented contig")
}

rm(FRAGMENTS)

#Extract chromosome info for assembly file
HAP12_assembly %>%
  filter(!grepl("^>",V1)) %>%
  mutate(CHR = row_number()) %>%
  gather(contig_order_in_chr, order_sign, V1:V100) %>%
  filter(!is.na(order_sign)) %>%
  mutate(contig_order_in_chr=as.numeric(gsub("V","",contig_order_in_chr))) %>%
  mutate(order = abs(as.numeric(order_sign))) %>%
  arrange(CHR,
          contig_order_in_chr)  %>%
  full_join(.,
            select(filter(HAP12_assembly,grepl("^>",V1)),fragment_ID = V1,original_order = V2),
            by=c("order"="original_order")) %>%
  mutate(fragment_ID = gsub(">","",fragment_ID)) %>%
  mutate(HAP = ifelse(grepl("h1",fragment_ID),"H1","H2")) -> mid_data

#to solve the problem of moving contigs between scaffolds of different haplotypes
mid_data %>%
  group_by(CHR,HAP) %>%
  tally() %>%
  ungroup() %>%
  group_by(CHR) %>%
  filter(n==max(n)) %>%
  ungroup() %>%
  select(CHR,HAP_new = HAP) -> CHRHAP

  mid_data %>%
  full_join(.,
            CHRHAP) %>%
  group_by(HAP_new) %>%
    arrange(HAP_new,CHR) %>%
    mutate(NEW_CHR = as.integer(factor(CHR))) %>%
  ungroup() %>%
  mutate(Chromosome = paste(HAP_new,"HiC_scaffold",NEW_CHR,sep = "_")) %>%
  mutate(orientation = ifelse(order_sign<0,"R","F")) %>%
  select(Chromosome,
         contig_order_in_chr,
         fragment_ID,
         orientation) -> CONTIG_ORDER_CHR

rm(mid_data,
CHRHAP)

HAP12_FASTA %>%
  full_join(.,
            CONTIG_ORDER_CHR,by=c("name" = "fragment_ID")) -> REV_FOR_CONTIG

rm(#CONTIG_ORDER_CHR,
   HAP12_FASTA,
   HAP12_assembly)
#only fragments that needs to be reverse complemented
REV_FOR_CONTIG %>%
  filter(orientation == "R") %>%
  select(-orientation) -> REV_CONTIG
#only fragments that does needs to be reverse complemented
REV_FOR_CONTIG %>%
  filter(orientation == "F") %>%
  select(-orientation) -> FOR_CONTIG
#reverse complement the fragments
REV_COMP_CONTIGS <- NULL
for (CONTIG in 1:nrow(REV_CONTIG)) {
  REV_CONTIG %>%
    filter(name == as.character(REV_CONTIG[CONTIG,2])) %>%
    pull(sq) %>%
    reverse() %>%
    complement() %>%
    enframe() %>%
    mutate(name = as.character(REV_CONTIG[CONTIG,2]),
           Chromosome = as.character(REV_CONTIG[CONTIG,3]),
           contig_order_in_chr = as.numeric(as.character(REV_CONTIG[CONTIG,4]))) %>%
    select(sq = value,
           name,
           Chromosome,
           contig_order_in_chr) %>%
    bind_rows(.,
              REV_COMP_CONTIGS) -> REV_COMP_CONTIGS
}

rm(REV_FOR_CONTIG,
   REV_CONTIG,
   CONTIG)

#merge the reversed complemented fragments with the rest corrected contig data set
bind_rows(FOR_CONTIG,
          REV_COMP_CONTIGS) -> CORRECTED_CONTIGS

rm(FOR_CONTIG,
   REV_COMP_CONTIGS)

CORRECTED_CONTIGS %>%
  distinct(Chromosome) %>%
  pull(Chromosome) -> CHROMOSOME


#assemble fragments to generate scaffolds and add 100Ns between contigs to build chromosomes
FINAL_FASTA <- NULL
for (CHR in 1:length(CHROMOSOME)) {
  CORRECTED_CONTIGS %>%
    filter(Chromosome == as.character(CHROMOSOME[CHR])) -> CHR_CONTIGS

  sq(rep(strrep("N",100),nrow(CHR_CONTIGS)),alphabet = "dna_ext") %>%  #100N GAP
    enframe() %>%
    mutate(contig_order_in_chr = seq(0.5,nrow(CHR_CONTIGS),1),
           name = "GAP",
           Chromosome = as.character(CHROMOSOME[CHR])) %>%
    select(sq = value,
           name,
           Chromosome,
           contig_order_in_chr) %>%
    bind_rows(.,
              CHR_CONTIGS) %>%
    arrange(contig_order_in_chr) %>% 
    filter(contig_order_in_chr != 0.5) %>%
    pull(sq) %>%
    collapse() %>%
    enframe() %>%
    mutate(name = as.character(CHROMOSOME[CHR])) %>%
    select(value,
           name) %>%
    bind_rows(.,
              FINAL_FASTA) -> FINAL_FASTA
  
  rm(CHR_CONTIGS)
}

rm(CHR,
   CHROMOSOME)

FINAL_FASTA %>%
  mutate(length = get_sq_lengths(value)) %>%
  arrange(desc(length)) %>%
  select(-length)-> FINAL_FASTA

#save assembled contigs for haplotype 1
write_fasta(pull(filter(FINAL_FASTA,grepl("^H1",name)),value),pull(filter(FINAL_FASTA,grepl("^H1",name)),name),paste0(SAVE_DIR,"/",PREFIX,"_hap1.reviewed.chr_assembled.fasta"))

#save assembled contigs for haplotype 2
write_fasta(pull(filter(FINAL_FASTA,grepl("^H2",name)),value),pull(filter(FINAL_FASTA,grepl("^H2",name)),name),paste0(SAVE_DIR,"/",PREFIX,"_hap2.reviewed.chr_assembled.fasta"))

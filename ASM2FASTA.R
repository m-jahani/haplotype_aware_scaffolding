library(tidyverse)
library(data.table)
library(tidysq)


####################################################################Read Data#####################################################################   
# read *reviewed assembly hap1_hap2 ()assembly file after juicebox curation
read.table("/Users/mojtabajahani/Downloads/AGA10_r0/mixedhap12.final.review.assembly",
      fill = TRUE, 
      col.names = paste("V", 1:100, sep = "")# with the assumption that each chromosome does not contain more than 100 conigs 
) -> HAP12_assembly


read_fasta("/Users/mojtabajahani/Downloads/AGA10_r0/AGA10.hic.hap12.p_ctg.fasta") -> HAP12_FASTA #merged hap1 and hap2 contigs
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
  arrange(as.numeric(Fragment_order))%>%
  ungroup() %>%
  mutate(end = ifelse(Fragment_order==1,length,lag(length)+length)) %>%
  mutate(start = (end-length)+1) %>% 
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

#remove the the original contig and replce the fragments
HAP12_FASTA %>% 
  filter(!name %in%  pull(distinct(coord_fragments,contig))) %>% 
  bind_rows(.,
            FRAGMENTED_CONTIGS) -> HAP12_FASTA_FRAG
}else {
  print("NO fragmented contig")
}

rm(FRAGMENTS,
   HAP12_FASTA,
   coord_fragments,
   FRAGMENTED_CONTIGS,
   CONTIG)


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
  mutate(HAP=ifelse(grepl("h1",fragment_ID),"H1","H2")) %>% 
  group_by(HAP) %>%
    arrange(HAP,CHR)  %>%
    mutate(NEW_CHR = as.integer(factor(CHR))) %>%
  ungroup() %>%
  mutate(Chromosome=paste(HAP,"HiC_scaffold",NEW_CHR,sep = "_")) %>%
  mutate(orientation=ifelse(order_sign<0,"R","F")) %>% 
  select(Chromosome,
         contig_order_in_chr,
         fragment_ID,
         orientation) -> CONTIG_ORDER_CHR 


HAP12_FASTA_FRAG %>% 
  full_join(.,
            CONTIG_ORDER_CHR,by=c("name" = "fragment_ID")) -> REV_FOR_CONTIG

rm(CONTIG_ORDER_CHR,
   HAP12_FASTA_FRAG,
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

#assemble fragments to generate scaffolds
FINAL_FASTA <- NULL
for (CHR in 1:length(CHROMOSOME)) {
  CORRECTED_CONTIGS %>%
    filter(Chromosome == as.character(CHROMOSOME[CHR])) %>% 
    arrange(contig_order_in_chr) %>%
    pull(sq) %>%
    collapse() %>%
    enframe() %>%
    mutate(name = as.character(CHROMOSOME[CHR])) %>%
    select(sq = value,
           name) %>%
    bind_rows(.,
              FINAL_FASTA) -> FINAL_FASTA
}

rm(CHR,
   CHROMOSOME)

FINAL_FASTA %>% 
  mutate(length = get_sq_lengths(sq)) %>% 
  arrange(desc(length)) %>%
  select(-length)-> FINAL_FASTA
  
#save assembled contigs
write_fasta(pull(FINAL_FASTA,sq),pull(FINAL_FASTA,name),"/Users/mojtabajahani/Downloads/AGA10_r0/test.fasta")

rm(FINAL_FASTA)


#save contig sequences
CORRECTED_CONTIGS %>% 
  select(-Chromosome,
         -contig_order_in_chr) %>%
  arrange(name) -> CORRECTED_CONTIGS_sq
  
write_fasta(pull(CORRECTED_CONTIGS_sq,sq),pull(CORRECTED_CONTIGS_sq,name),"/Users/mojtabajahani/Downloads/AGA10_r0/test_contig.fasta")

rm(CORRECTED_CONTIGS_sq)
#position of contigs in each chromosome
CORRECTED_CONTIGS %>% 
  mutate(length = get_sq_lengths(sq)) %>% 
  select(-sq) %>%
  group_by(Chromosome) %>%
  arrange(Chromosome,contig_order_in_chr) %>% 
  mutate(end = cumsum(length) ) %>% 
  ungroup() %>% 
  mutate(start = ifelse(contig_order_in_chr == 1,1,lag(end))) %>%
  select(Chromosome,
         start,
         end,
         contig = name) %>%
  fwrite("/Users/mojtabajahani/Downloads/AGA10_r0/test",
         sep = "\t",
         quote = F,
         col.names = F)

rm(CORRECTED_CONTIGS)
  
  
  
library(tidyverse)
library(tidysq)
library(data.table)

options(scipen = 9999)
args = commandArgs(trailingOnly = TRUE)

# MIX_ASSEM <- args[1]
# bed_file <- args[2]
# PREFIX <- args[3]

MIX_ASSEM <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/MANUALCURATION/HiC_map/AGA10.hic.hap1.p_ctg_AGA10.hic.hap2.p_ctg.review.assembly"
bed_file <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/MANUALCURATION/PLOT/edit1_asm.csv" #two culmns, scaffold name:(H1_HiC_scaffold_10), location to cut on scaffold: 53000
PREFIX <- "Edit1"
####################################################################Read Data#####################################################################   
# read *reviewed assembly hap1_hap2 ()assembly file after juicebox curation
read.table(MIX_ASSEM,
      fill = TRUE, 
      col.names = paste("V", 1:100, sep = "")# with the assumption that each chromosome does not contain more than 100 conigs 
) -> HAP12_assembly

fread(bed_file) %>%
  arrange(V1,as.numeric(V2)) -> bed
##################################################Contigs/fragment/Debris order in assembly file####################################################          
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
            select(filter(HAP12_assembly,grepl("^>",V1)),fragment_ID = V1,original_order = V2, length = V3),
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
  group_by(Chromosome) %>%
  mutate(end = cumsum(length)) %>%
  mutate(start = lag(end)) %>%
  mutate(start = replace_na(start,1)) %>%
  ungroup() %>% 
  separate(fragment_ID, 
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>%
  separate(lable1,
           into = c("Fragment","Fragment_order")) %>% 
  select(Chromosome,
         fragment_ID,
         contig,
         Fragment,
         Fragment_order,
         contig_order_in_chr,
         orientation,
         start,
         end,
         length,
         order) -> CONTIG_ORDER_CHR 


for (i in 1:nrow(bed)) { #1:nrow(bed)
  
CONTIG_ORDER_CHR %>%
  filter(Chromosome == as.character(bed[i,1])) %>%
  filter(start <= as.numeric(bed[i,2])) %>%
  filter(end >= as.numeric(bed[i,2]))  -> target_contig

target_contig %>%
  mutate(Fragment_order = as.numeric(replace_na(Fragment_order,"1"))) %>%
  mutate(Fragment = replace_na(Fragment,"fragment")) %>%
  mutate(end = as.numeric(bed[i,2]) ) %>% 
  mutate(length = ifelse(as.numeric(start)==1,
                         (as.numeric(end) - as.numeric(start))+1,
                         (as.numeric(end) - as.numeric(start)))) -> row_one

target_contig %>%
  mutate(Fragment = "fragment") %>%
  mutate(Fragment_order = ifelse(is.na(Fragment_order),
    as.numeric(row_one[1,5])+1,
   max(
    pull(
      mutate(
        mutate(
          filter(
            CONTIG_ORDER_CHR,
            contig == as.character(target_contig[1,3])),
          Fragment_order = replace_na(Fragment_order,"0")),
        Fragment_order = as.numeric(Fragment_order)),
      Fragment_order))+1)) %>%
  mutate(start = as.numeric(bed[i,2])) %>%
  mutate(contig_order_in_chr = as.numeric(row_one[1,6])+1) %>%
  mutate(order = as.numeric(row_one[1,11])+1) %>%
  mutate(length = (as.numeric(end) - as.numeric(start))) -> row_two

rbind(row_one,
      row_two) %>%
  mutate(fragment_ID = paste0(contig,":::",Fragment,"_",Fragment_order)) -> edited_fragments

sum(edited_fragments$length)

CONTIG_ORDER_CHR %>%
  filter(Chromosome == as.character(bed[i,1])) %>%
  filter(fragment_ID != as.character(target_contig[1,2])) %>%
  mutate(contig_order_in_chr = as.numeric(contig_order_in_chr))  %>%
  mutate(contig_order_in_chr = ifelse(contig_order_in_chr < as.numeric(max(edited_fragments[,6])) , contig_order_in_chr,contig_order_in_chr+1)) %>%
  arrange(contig_order_in_chr)  %>%
  rbind(.,
  filter(CONTIG_ORDER_CHR,Chromosome != as.character(bed[i,1]))) %>%
  mutate(order = as.numeric(order))  %>%
  mutate(order = ifelse(order < as.numeric(max(edited_fragments[,11])) , order,order+1)) %>%
  rbind(.,
        edited_fragments) %>%
  separate(Chromosome,into = c("HAP","num"),sep = "_HiC_scaffold_",remove=F) %>%
  arrange(as.numeric(num),HAP,contig_order_in_chr) %>%
  select(-HAP,-num)  ->  CONTIG_ORDER_CHR

rm(target_contig,row_one,row_two,edited_fragments)
}

CONTIG_ORDER_CHR %>% 
  arrange(order) %>%
  mutate(V1=paste0(">",fragment_ID)) %>% 
  select(V1,
         order,
         length) -> HAP1_HAP2_CONTOGS_NEW_ORDER
  
CONTIG_ORDER_CHR %>% 
  mutate(order_sign=ifelse(orientation == "F",order, -1 * as.numeric(order))) %>%
  select(Chromosome,
         contig_order_in_chr,
         order_sign) %>%
  mutate(order_sign=as.character(order_sign)) %>%
  spread(contig_order_in_chr,order_sign) %>%
  separate(Chromosome,into = c("HAP","num"),sep = "_HiC_scaffold_",remove = T) %>%
  arrange(as.numeric(num),HAP) %>%
  select(-num,-HAP) -> CHR_ORDER_HAP12

  
bind_rows(mutate_if(setNames(HAP1_HAP2_CONTOGS_NEW_ORDER,c(1,2,3)),is.numeric, as.character),CHR_ORDER_HAP12)  %>%
  fwrite(paste0(gsub("review.assembly","",MIX_ASSEM),
                PREFIX,
                ".review.assembly"),
         col.names = F,
         sep = " ")
  



  
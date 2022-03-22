library(data.table)
library(tidyverse)


args = commandArgs(trailingOnly = TRUE)

AASEM1 <- args[1] #reviewed assembly hap1
AASEM2 <- args[2] #reviewed assembly hap2
map_hap12 <- args[3] #read minimap2 alignment between hap1 and hap2 contigs in paf format
####################################################################Read Data#####################################################################   
# read *reviewed assembly hap1
read.table(AASEM1,
           fill = TRUE, 
           col.names = paste("V", 1:100, sep = "")# with the assumption that each chromosome does not contain more than 100 conigs 
           ) -> HAP1

# read *reviewed assembly hap2
read.table(AASEM2, 
           fill = TRUE, 
           col.names = paste("V", 1:100, sep = "")) -> HAP2

# read minimap2 alignment between hap1 and hap2 contigs
read.table(map_hap12,
                       fill = TRUE) -> HAP1_HAP2_ALIGNMENT
##################################################Contigs/fragment/Debris order in assembly file####################################################          

#calculate the original length of fragmented contigs in Hap1 (for sorting) 
HAP1 %>%
  filter(grepl("^>",V1)) %>% 
  select(fragment_ID = V1,
         original_order = V2,
         length = V3)  %>%
  separate(fragment_ID, 
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>% 
  group_by(contig) %>% 
  summarize(total_length=sum(length)) %>%
  ungroup() -> CONTIG_LENGTH_HAP1

#prepare contig data for hap 1    
HAP1 %>%
  filter(grepl("^>",V1)) %>% 
  select(fragment_ID = V1,
         original_order = V2,
         length = V3)  %>%
  separate(fragment_ID, 
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>% 
  full_join(.,
            CONTIG_LENGTH_HAP1) %>%
  mutate(HAP = 1) %>%
  select(fragment_ID,
         original_order,
         length,
         total_length,
         HAP) -> DATA_HAP1

rm(CONTIG_LENGTH_HAP1)

#calculat the originl length of fragmented contigs in Hap2 (for sorting) 
HAP2 %>%
  filter(grepl("^>",V1)) %>% 
  select(fragment_ID = V1,
         original_order = V2,
         length = V3)  %>%
  separate(fragment_ID, 
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>%
  group_by(contig) %>%
  summarize(total_length=sum(length)) %>%
  ungroup() -> CONTIG_LENGTH_HAP2

#prepare contig data for hap 2
HAP2 %>%
  filter(grepl("^>",V1)) %>% 
  select(fragment_ID = V1,
         original_order = V2,
         length = V3)  %>%
  separate(fragment_ID, 
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>% 
  full_join(.,
            CONTIG_LENGTH_HAP2) %>%
  mutate(HAP = 2) %>%
  select(fragment_ID,
         original_order,
         length,
         total_length,
         HAP) -> DATA_HAP2

rm(CONTIG_LENGTH_HAP2)

#hap1 and hap 2 contigs with the new order
rbind(DATA_HAP1,
      DATA_HAP2) %>%
  arrange(
    desc(total_length)) %>% 
  mutate(new_order = row_number()) %>%
  select(fragment_ID,
         new_order,
         length) -> HAP1_HAP2_CONTOGS_NEW_ORDER

##################################################Contig orders with mixed haps in each chromosome####################################################

#new order to original order
rbind(DATA_HAP1,
      DATA_HAP2) %>%
  arrange(
    desc(total_length)) %>%
  mutate(new_order = row_number()) %>% 
  select(HAP,
         original_order,
         new_order) -> ORDER

#Extract chromosome info for hap1
HAP1 %>% 
  filter(!grepl("^>",V1)) %>%
  mutate(CHR = row_number()) %>%
  gather(contig_order_in_chr, original_order_sign, V1:V100) %>%
  filter(!is.na(original_order_sign)) %>%
  mutate(contig_order_in_chr=as.numeric(gsub("V","",contig_order_in_chr))) %>%
  mutate(original_order = abs(as.numeric(original_order_sign))) %>%
  arrange(CHR,
          contig_order_in_chr) %>% 
  mutate(HAP = 1) -> CHR_CONTIG_HAP1_INFO

#Calculate chromosome length for hap1
CHR_CONTIG_HAP1_INFO %>% 
  full_join(.,
            select(
              filter(
                HAP1,grepl("^>",V1)),
              original_order = V2,
              length = V3) ) %>% 
  group_by(HAP,CHR) %>%
  summarise(CHR_length=sum(length)) %>% 
  ungroup()  -> CHR_LENGTH_HAP1

#Extract chromosome info for hap2
HAP2 %>%
  filter(!grepl("^>",V1)) %>%
  mutate( CHR = row_number()) %>%
  gather(contig_order_in_chr, original_order_sign, V1:V100) %>%
  filter(!is.na(original_order_sign)) %>%
  mutate(contig_order_in_chr=as.numeric(gsub("V","",contig_order_in_chr))) %>%
  mutate(original_order=abs(as.numeric(original_order_sign))) %>%
  arrange(CHR,
          contig_order_in_chr) %>% 
  mutate(HAP = 2) -> CHR_CONTIG_HAP2_INFO

#Calculate chromosome length for hap2
CHR_CONTIG_HAP2_INFO %>%
  full_join(.,
            select(
              filter(
                HAP2,grepl("^>",V1)),
              original_order = V2,
              length = V3) ) %>% 
  group_by(HAP,CHR) %>%
  summarise(CHR_length=sum(length)) %>% 
  ungroup()  -> CHR_LENGTH_HAP2

#Merge chromosome length for hap1 and hap2
rbind(CHR_LENGTH_HAP1,CHR_LENGTH_HAP2) -> CHR_LENGTH_HAP12

rm(CHR_LENGTH_HAP1,
   CHR_LENGTH_HAP2,
   HAP1,
   HAP2)


#######################################################Chromosome order in contact map#########################################################
#corresponding contigs in hap1 
HAP1_HAP2_ALIGNMENT %>%
  mutate(alignment_perc = V11/V2) %>%
  filter(V12 > 1) %>%
  group_by(V1) %>%
  filter(alignment_perc == max(alignment_perc)) %>%
  ungroup() %>%
  distinct(V1,
           V6) %>%
  mutate(contig = paste0(">",V1),
         corresponding = paste0(">",V6)) %>%
  select(contig,
         corresponding) -> REF_ASM

#corresponding contigs in hap2
HAP1_HAP2_ALIGNMENT %>%
  mutate(alignment_perc = V11/V7) %>%
  filter(V12 > 1) %>%
  group_by(V6) %>%
  filter(alignment_perc == max(alignment_perc)) %>%
  ungroup() %>%
  distinct(V1,
           V6) %>%
  mutate(contig = paste0(">",V6),
         corresponding = paste0(">",V1)) %>%
  select(contig,
         corresponding) -> ASM_REF

rbind(REF_ASM,ASM_REF) -> H12_CORES_CONTIGS

rm(REF_ASM,ASM_REF,HAP1_HAP2_ALIGNMENT)

#match corresponding chromosomes in haps
rbind(DATA_HAP1,
      DATA_HAP2) %>%
  arrange(
    desc(total_length)) %>%
  mutate(new_order = row_number()) %>%
  separate(fragment_ID,
           into = c("contig"),
           sep = ":::",
           remove=F) %>%
  full_join(.,H12_CORES_CONTIGS) -> CONTIG_DATA

rm(DATA_HAP1,
   DATA_HAP2,
   H12_CORES_CONTIGS)

CONTIG_DATA %>%
  select(HAP,
         new_order,
         corresponding) %>%
  full_join(.,
            select(CONTIG_DATA,
                   corresponding = contig,
                   new_order_cor = new_order)
            ) %>%
  select(-corresponding) %>%
  filter(!is.na(new_order_cor),
         !is.na(HAP)) -> CONTIG_COR

rm(CONTIG_DATA)

#order_id of chromosome pairs
rbind(CHR_CONTIG_HAP1_INFO,
      CHR_CONTIG_HAP2_INFO)  %>%
  select(CHR,
         original_order,
         HAP) %>%
  full_join(.,
            ORDER) %>%
  select(-original_order) %>%
  full_join(.,
            CONTIG_COR) -> DATA_CHR_PAIR

rm(CONTIG_COR)

DATA_CHR_PAIR %>% filter(HAP==1) %>% select(CHR_HAP1 = CHR, new_order_hap1 = new_order , new_order_hap2 = new_order_cor) -> DATA_CHR_PAIR_HAP1
DATA_CHR_PAIR %>% filter(HAP==2) %>% select(CHR_HAP2 = CHR, new_order_hap2 = new_order , new_order_hap1 = new_order_cor) -> DATA_CHR_PAIR_HAP2

rm(DATA_CHR_PAIR)

 full_join(DATA_CHR_PAIR_HAP1,DATA_CHR_PAIR_HAP2) %>%
   group_by(CHR_HAP1,CHR_HAP2) %>%
   tally() %>%
   ungroup() %>%
   filter(!is.na(CHR_HAP1),
          !is.na(CHR_HAP2)) %>%
   arrange(desc(n)) %>%
   left_join(.,
             select(filter(CHR_LENGTH_HAP12,HAP == 1), CHR_HAP1 = CHR ,CHR_length_HAP1 = CHR_length)) %>%
   left_join(.,
             select(filter(CHR_LENGTH_HAP12,HAP == 2), CHR_HAP2 = CHR ,CHR_length_HAP2 = CHR_length)) %>%
   mutate(length_diff_kb = abs(CHR_length_HAP2-CHR_length_HAP1)/1000) %>%
   filter(CHR_length_HAP1 > 10000000, CHR_length_HAP2 > 10000000) %>%
   group_by(CHR_HAP1) %>%
   filter(length_diff_kb == min(length_diff_kb)) %>%
   ungroup() %>%
   group_by(CHR_HAP2) %>%
   filter(length_diff_kb == min(length_diff_kb)) %>%
   ungroup() %>%
   group_by(CHR_HAP1) %>%
   filter(n == min(n)) %>%
   ungroup() %>%
   group_by(CHR_HAP2) %>%
   filter(n == min(n)) %>%
   ungroup() %>%
   arrange(CHR_length_HAP1)  %>%
   select(CHR_HAP1,
          CHR_HAP2) %>%
   mutate(chr_ord = row_number()) %>%
   gather(HAP,CHR,CHR_HAP1,CHR_HAP2) %>%
   mutate(HAP=as.numeric(gsub("CHR_HAP","",HAP))) %>%
   arrange(chr_ord,HAP) %>%
   mutate(chr_hap_ord = row_number()) %>%
   select(CHR,
          HAP,
          chr_hap_ord) -> CHR_HAP_ORDER
 
 rm(DATA_CHR_PAIR_HAP1,
    DATA_CHR_PAIR_HAP2)
################################################################Building the assembly file for 3D-DNA######################################################

#re-build the chr orders and length with merged haplotype
rbind(CHR_CONTIG_HAP1_INFO,
      CHR_CONTIG_HAP2_INFO)  %>% 
  full_join(.,
            CHR_LENGTH_HAP12) %>% 
  full_join(.,
            ORDER) %>% 
  mutate(new_order_sign=if_else(as.numeric(original_order_sign) > 0,
                                (new_order*1),
                                (new_order*-1)
  )
  )  %>% 
  select(HAP,
         CHR,
         CHR_length,
         contig_order_in_chr,
         new_order_sign) %>% 
  spread(contig_order_in_chr,new_order_sign) %>%
  left_join(.,CHR_HAP_ORDER) %>%
  arrange(chr_hap_ord,
          desc(CHR_length))  %>%
   select(-HAP,
          -CHR,
          -CHR_length,
          -chr_hap_ord) %>%
   mutate_if(is.numeric, as.character) -> CHR_ORDER_HAP12
 
 rm(ORDER,
    CHR_CONTIG_HAP1_INFO,
    CHR_CONTIG_HAP2_INFO,
    CHR_LENGTH_HAP12,
    CHR_HAP_ORDER)
 
 bind_rows(mutate_if(setNames(HAP1_HAP2_CONTOGS_NEW_ORDER,c(1,2,3)),is.numeric, as.character),CHR_ORDER_HAP12)  %>%
   fwrite("/Users/mojtabajahani/Downloads/AGA10_r0/mixedhap12.final.review.assembly",
          col.names = F,
          sep = " ")

 
rm(HAP1_HAP2_CONTOGS_NEW_ORDER,
   CHR_ORDER_HAP12)


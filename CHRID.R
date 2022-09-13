library(tidyverse)
library(tidysq)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

MAP <- args[1] #minimap result in PAF format
GAP_FASTA <- args[2]
PREFIX <- args[3] #AGA10_Hap1
SAVE_DIR <- args[4]

GAP_FASTA <- "/DATA/home/mjahani/CURATION/AGA10/ASM2FASTA/FINAL_GAP/AGA10_FINAL_GAP_hap1.reviewed.chr_assembled.fasta"
MAP <- "/DATA/home/mjahani/CURATION/AGA10/minimap/AGA10_FINAL_GAP_hap1.reviewed.chr_assembled_GCF_900626175.2_cs10_10CHR.fa.paf" 
PREFIX <- "AGA10_Hap1"
SAVE_DIR <- "/DATA/home/mjahani/CURATION"



read.table(MAP,fill = TRUE) %>%
  group_by(V1,V6) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(V1) %>% 
  filter(n==max(n)) %>% 
  ungroup() %>%
  select(CS10 = V1,
        name = V6) %>%
  full_join(.,
            data.frame(CS10 = c("NC_044371.1",
                              "NC_044375.1",
                              "NC_044372.1",
                              "NC_044373.1",
                              "NC_044374.1",
                              "NC_044377.1",
                              "NC_044378.1",
                              "NC_044379.1",
                              "NC_044376.1",
                              "NC_044370.1"),
                       CHROMOSOME = c("Chr01",
                                      "Chr02",
                                      "Chr03",
                                      "Chr04",
                                      "Chr05",
                                      "Chr06",
                                      "Chr07",
                                      "Chr08",
                                      "Chr09",
                                      "ChrX"))) %>%
  mutate(CHROMOSOME = paste0(PREFIX,"_",CHROMOSOME)) %>%
  select(name,
         CHROMOSOME) -> Conversion_table


read_fasta(GAP_FASTA) %>%
  full_join(.,
            Conversion_table) %>%
  mutate(row=row_number()-10) %>%
  mutate(CHROMOSOME = ifelse(is.na(CHROMOSOME),
                             paste0(PREFIX,"_","Contig",row),
                             CHROMOSOME)) %>%
  arrange(CHROMOSOME) %>%
  mutate(order = as.numeric(ifelse(as.numeric(row)<1,row_number()-10,row))) %>%
  arrange(order)  %>%
  select(sq,
         CHROMOSOME) -> FINAL_FASTA

write_fasta(pull(FINAL_FASTA,sq),pull(FINAL_FASTA,CHROMOSOME),paste0(SAVE_DIR,"/",PREFIX,".GAP.CHR_ID.reviewed.chr_assembled.fasta"))


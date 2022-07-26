library(data.table)
library(tidyverse)

options(scipen = 9999)

args = commandArgs(trailingOnly = TRUE)

coordinates_file <- args[1]
CONTIG = args[2]

fread(coordinates_file) -> data

data %>% 
  group_by(V1) %>%
  summarise(scaf_length=max(V3)) %>% 
  separate(V1, 
           into = c("hap","scaf_num"), 
           sep="_HiC_scaffold_", 
           remove = F) %>%
  arrange(as.numeric(scaf_num)) %>%
  mutate(length = lag(scaf_length)) %>%
  mutate(length = replace_na(length, 0)) %>%
  mutate(cum_length = cumsum(length)) %>%
  select(V1,
         cum_length) -> scf_length
  

print(data %>% 
  full_join(.,scf_length) %>% 
  mutate(map_start = as.numeric(V2) + as.numeric(cum_length)) %>%
  mutate(map_end = as.numeric(V3) + as.numeric(cum_length)) %>% 
  select(scaffold = V1,
         contig = V4,
         map_start,
         map_end) %>% 
  filter(contig == CONTIG))




  
  

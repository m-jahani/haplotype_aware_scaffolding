library(data.table)
library(tidyverse)


args = commandArgs(trailingOnly = TRUE)
THREE_PRIM <- args[1]
FIVE_PRIM <- args[2]
LINKAGE_MAP <- args[3]
SAVE_DIR <- args[4]

fread(THREE_PRIM) %>% 
  select(QNAME = V1,
         RNAME_3 = V2,
         POS_3 = V3) %>% 
  separate(QNAME, into = c("Marker","extra"), sep = "::", remove = T) %>% 
  select(-extra) -> three

fread(FIVE_PRIM) %>% 
  select(QNAME = V1,
         RNAME_5 = V2,
         POS_5 = V3) %>% 
  separate(QNAME, into = c("Marker","extra"), sep = "::", remove = T) %>% 
  select(-extra) %>% 
  full_join(.,three)  %>% 
  filter(!is.na(POS_3)) %>% #filter out 3prime markers without corresponding 5prime marker
  filter(!is.na(POS_5)) %>% #filter out 5prime markers without corresponding 3prime marker
  mutate(same=RNAME_5==RNAME_3) %>% 
  filter(same != FALSE) %>% #filter out markers that each boundary mapped on a different scaffold
  select(Marker,scaffold = RNAME_5,POS_5,POS_3) -> mapped_markers

fread(LINKAGE_MAP) %>%
  inner_join(.,mapped_markers) %>%
  mutate(ID=paste0(`Linkage group`,"::",`Genetic distance (cM)`)) %>% 
  mutate(BP_position = rowMeans(select(., starts_with("POS")), na.rm = TRUE)) %>%  #mean value of position 3 prime and 5 prime as marker position
  select(Marker,
         scaffold,
         BP_position,
         linkage_group = `Linkage group`,
         cM = `Genetic distance (cM)`) -> genetic_map_AK_HAP2

genetic_map_AK_HAP2 %>%
  group_by(scaffold) %>% 
  arrange(BP_position,.by_group = TRUE) %>%
  mutate(delta_cM=cM-lag(cM, default = first(cM)),
         delta_bp=BP_position-lag(BP_position, default = first(BP_position)))  %>% 
  ungroup() %>%
  mutate(recomb_rate=abs(delta_cM)/abs(delta_bp)) %>% 
  select(scaffold,BP_position,recomb_rate)  %>% 
  filter(!is.na(recomb_rate)) %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                gsub(".sam","",gsub("3_primeBoundary_","",gsub(".*/","",THREE_PRIM))),
                ".recombination"),
         sep = "\t",
         col.names = F,
         quote = F)

genetic_map_AK_HAP2 %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                gsub(".sam","",
                     gsub("3_primeBoundary_","",
                          gsub(".*/","",THREE_PRIM))),
                ".GeneticMap"),
         sep = "\t",
         col.names = F,
         quote = F)
  




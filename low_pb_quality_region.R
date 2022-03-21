library(data.table)
library(tidyverse)
library(IRanges)

args = commandArgs(trailingOnly = TRUE)

BED <- args[1] #low quality regions from Merqury

fread(BED) -> LOW_Q_BED

LOW_Q_BED %>% distinct(V1) %>% pull(V1) -> SCAFF
MERGED_LOW_Q <- NULL
for (i in 1:length(SCAFF)) {
  LOW_Q_BED %>% 
    filter(V1 == SCAFF[i]) %>%
    select(-V1) -> coordinates
  
    IRanges::reduce(
      IRanges(pull(coordinates,V2),pull(coordinates,V3))) %>% 
      as.data.frame() %>%
      mutate(scaffold = SCAFF[i]) %>%
      select(scaffold,
             start,
             end) %>%
      rbind(.,MERGED_LOW_Q) -> MERGED_LOW_Q
}

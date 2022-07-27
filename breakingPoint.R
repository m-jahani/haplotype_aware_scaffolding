library(data.table)
library(tidyverse)
library(pafr)

options(scipen = 9999)
args = commandArgs(trailingOnly = TRUE)

Synteny <- args[1]
save_dir <- args[2]

read_paf(Synteny) %>%
  as.data.frame() %>% 
  filter(mapq > 0) %>% 
  mutate(new_strand=ifelse(strand=="+","Forward","Reverse"))%>%
  mutate(new_tstart= ifelse(new_strand == "Forward",tstart,tend)) %>% 
  mutate(new_tend= ifelse(new_strand == "Forward",tend,tstart)) %>% 
  select(qname,
         qstart,
         qend,
         tname,
         tstart = new_tstart,
         tend = new_tend,
         new_strand) %>%
  separate(qname,
           into = c("hap","qchr_num"),
           sep = "_HiC_scaffold_",
           remove=F) %>%
  filter(as.numeric(qchr_num) %in% seq(1,10)) %>%
  select(-hap) %>%
  separate(tname,
           into = c("hap","tchr_num"),
           sep = "_HiC_scaffold_",
           remove=F) %>%
  filter(as.numeric(tchr_num) %in% seq(1,10)) %>%
  mutate(test=as.numeric(qchr_num)  == as.numeric(tchr_num)) %>% 
  filter(test == "TRUE") %>% 
  select(-tchr_num,
         -qchr_num,
         -test) %>%
  mutate(new_qstart=ifelse(as.numeric(qstart) < as.numeric(qend),qstart,qend)) %>%
  mutate(new_qend=ifelse(as.numeric(qstart) < as.numeric(qend),qend,qstart)) %>%
  mutate(new_tstart=ifelse(as.numeric(tstart) < as.numeric(tend),tstart,tend)) %>%
  mutate(new_tend=ifelse(as.numeric(tstart) < as.numeric(tend),tend,tstart)) %>%
  select(-qstart,-qend,-tstart,-tend) %>%
  select(qname,
         qstart = new_qstart,
         qend = new_qend,
         tname,
         tstart = new_tstart,
         tend = new_tend) %>%
  fwrite(paste0(save_dir,
                "/",
                gsub("paf","alignment.csv",
                     gsub(".*/","",Synteny))),
         col.names = T,
         sep = ",")




  
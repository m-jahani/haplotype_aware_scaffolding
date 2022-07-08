library(plotly)
library(data.table)
library(tidyverse)
library(pafr)
options(scipen = 9999)
scaffold_H1 <- "H1_HiC_scaffold_1"
scaffold_H2 <- "H2_HiC_scaffold_1" 
WIND_SIZE=200000 #EDTA window size
##########################################################
coordinates_file <-"/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/AGA10hap12.reviwed_contig_chr_coord"
geneticmap_hap1 <-"/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/LINKAGE_MAP/AGA10hap1.reviewed.chr_assembled.GeneticMap"
geneticmap_hap2 <-"/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/LINKAGE_MAP/AGA10hap2.reviewed.chr_assembled.GeneticMap"
recomnination_hap1 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/LINKAGE_MAP/AGA10hap1.reviewed.chr_assembled.recombination"
recomnination_hap2 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/LINKAGE_MAP/AGA10hap2.reviewed.chr_assembled.recombination"
telomere_hap1 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/TELOMERE/AGA10hap1.reviewed.chr_assembled_telomeric_repeat_200000windows.bed"
telomere_hap2 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/TELOMERE/AGA10hap2.reviewed.chr_assembled_telomeric_repeat_200000windows.bed"
EDTA_hap1 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/EDTA/AGA10hap1.reviewed.chr_assembled_10chr.bed"
EDTA_hap2 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/EDTA/AGA10hap2.reviewed.chr_assembled_10chr.bed"
#DEPTH_hap1 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/depth/AGA10_hifi_hap1.sorted.smoothed.depth.bed"
#DEPTH_hap2 <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/depth/AGA10_hifi_hap2.sorted.smoothed.depth.bed"
original_coordinatess_for_EDTA <-"/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/AGA10hap12.reviwed_contig_chr_coord"
Synteny <- "/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/scaffold_result/SYNTENY/AGA10hap1.reviewed.chr_assembled_AGA10hap2.reviewed.chr_assembled.paf"

##########################################################Functions
#calculate the start and end position of fragmented contigs and their coordinates in scaffolds
fragmented_coordinates_modification <- function(ARGUMENT1) {
fread(ARGUMENT1) %>%  
  select(Scaffold = V1,
         contig = V4,
         start = V2,
         end = V3) %>% 
  separate(contig, 
           into = c("contig_base","fragment_ID","debris"), 
           sep = ":::", 
           remove = FALSE) -> contig_coor

contig_coor %>% 
  filter(!is.na(fragment_ID)) -> FRAGMENTS

if (nrow(FRAGMENTS) > 0 ) {
  
  FRAGMENTS %>%
    mutate(fragment_ID = as.numeric(gsub("fragment_","",fragment_ID))) %>% 
    group_by(contig_base) %>%
    arrange(fragment_ID) %>% 
    mutate(length = (as.numeric(end)-as.numeric(start))+1) %>% 
    mutate(end_contig = cumsum(length)) %>%
    mutate(start_contig = ifelse(fragment_ID == 1, 1, lag(end_contig)+1)) %>% 
    ungroup() %>%
    select(Scaffold,
           contig, #contig:::fargment
           contig_base, #contig base name (without fragment ID)
           start_scaffold = start,
           end_scaffold = end,
           start_contig,
           end_contig) -> FRAGMENTED_SE
  
  contig_coor %>% 
    filter(is.na(fragment_ID)) %>%
    mutate(end_contig = (as.numeric(end)-as.numeric(start))+1) %>%
    mutate(start_contig = 1) %>%
    select(Scaffold,
           contig,
           contig_base,
           start_scaffold = start,
           end_scaffold = end,
           start_contig,
           end_contig) %>%
    rbind(.,FRAGMENTED_SE) ->coordinates_internal
  
  
}else {
  
  contig_coor %>% 
    filter(is.na(fragment_ID)) %>%
    mutate(end_contig = (as.numeric(end)-as.numeric(start))+1) %>%
    mutate(start_contig = 1) %>%
    select(Scaffold,
           contig,
           contig_base,
           start_scaffold = start,
           end_scaffold = end,
           start_contig,
           end_contig) -> coordinates_internal
}

return(coordinates_internal)

rm(FRAGMENTS,
   FRAGMENTED_SE,
   contig_coor)
}

find the new position after reorder in coordinations-scaffolds to original contigs
scaffold2contig_ranges1 <- function(ARGUMENT1, #data (telomere, EDTA and etc)
                                   ARGUMENT2) { #original coordinates
  full_join(ARGUMENT1,
            ARGUMENT2) %>%
    mutate(start_overlap =
             ifelse(start >= start_scaffold  & start <= end_scaffold,
                    "YES",
                    "NO")) %>%
    mutate(end_overlap =
             ifelse(end >= start_scaffold  & end <= end_scaffold,
                    "YES",
                    "NO")) %>%
    mutate(start_end_overlap = ifelse(start_overlap == "NO" & end_overlap == "NO",
                                      "NO_OVERLAP",
                                      ifelse(start_overlap == "YES" & end_overlap == "NO",
                                             "START_OVERLAP",
                                             ifelse(start_overlap == "NO" & end_overlap == "YES",
                                                    "END_OVERLAP",
                                                    ifelse(start_overlap == "YES" & end_overlap == "YES",
                                                           "BOTH_OVERLAP",
                                                           "ERROR"))))) %>%
    select(-start_overlap,
           -end_overlap) %>%
    filter(start_end_overlap != "NO_OVERLAP",
           start_end_overlap != "ERROR",
           !is.na(start_end_overlap))  %>%
    mutate(start = ifelse(start_end_overlap == "END_OVERLAP", start_scaffold, start)) %>%
    mutate(end = ifelse(start_end_overlap == "START_OVERLAP", end_scaffold, end)) %>%
    mutate(start_contig_new = (((as.numeric(start) - as.numeric(start_scaffold))+1)+ as.numeric(start_contig))-1 ) %>%
    mutate(end_contig_new = (((as.numeric(end) - as.numeric(start_scaffold))+1)+ as.numeric(start_contig))-1 ) %>%
    select(contig_base = contig_base,
           start = start_contig_new,
           end  = end_contig_new,
           value) -> internal_result
  return(internal_result)
}

#find the new position after reorder in coordinations-scaffolds to original contigs
scaffold2contig_ranges <- function(ARGUMENT1, #data (telomere, EDTA and etc)
                                   ARGUMENT2) { #original coordinates
  full_join(ARGUMENT1, 
            ARGUMENT2) %>% 
    mutate(start_overlap = 
             ifelse(start >= start_scaffold  & start <= end_scaffold,
                    "YES",
                    "NO")) %>%
    mutate(end_overlap = 
             ifelse(end >= start_scaffold  & end <= end_scaffold,
                    "YES",
                    "NO")) %>% 
    mutate(between = 
             ifelse(start < start_scaffold  & end > end_scaffold,
                    "YES",
                    "NO")) %>% 
    mutate(start_end_overlap = ifelse(start_overlap == "NO" & end_overlap == "NO" & between == "NO",
                                      "NO_OVERLAP",
                                      ifelse(start_overlap == "YES" & end_overlap == "NO",
                                             "START_OVERLAP",
                                             ifelse(start_overlap == "NO" & end_overlap == "YES",
                                                    "END_OVERLAP",
                                                    ifelse(start_overlap == "YES" & end_overlap == "YES",
                                                           "BOTH_OVERLAP",
                                                           ifelse(between == "YES",
                                                                  "BETWEEB_OVERLAP",
                                                                  "ERROR")))))) %>% 
    select(-start_overlap,
           -end_overlap) %>%
    filter(start_end_overlap != "NO_OVERLAP",
           start_end_overlap != "ERROR",
           !is.na(start_end_overlap))  %>% 
    mutate(start = ifelse(start_end_overlap == "END_OVERLAP" | start_end_overlap == "BETWEEB_OVERLAP", 
                          start_scaffold, 
                          start)) %>%
    mutate(end = ifelse(start_end_overlap == "START_OVERLAP" | start_end_overlap == "BETWEEB_OVERLAP", 
                        end_scaffold, 
                        end)) %>% 
    mutate(start_contig_new = (((as.numeric(start) - as.numeric(start_scaffold))+1)+ as.numeric(start_contig))-1 ) %>%
    mutate(end_contig_new = (((as.numeric(end) - as.numeric(start_scaffold))+1)+ as.numeric(start_contig))-1 ) %>%
    select(contig_base = contig_base,
           start = start_contig_new,
           end  = end_contig_new,
           value) -> internal_result
  return(internal_result)
}


#original contigs to new fragmented contig scaffold position
contig2scaffold_ranges <- function(ARGUMENT1,#result of scaffold2contig
                                   ARGUMENT2) { #new coordinates
  ARGUMENT1 %>%
    select(contig_base,
           start,
           end,
           value) %>%
    full_join(.,ARGUMENT2) %>% 
    mutate(start_overlap = 
             ifelse(start >= start_contig  & start <= end_contig,
                    "YES",
                    "NO")) %>%
    mutate(end_overlap = 
             ifelse(end >= start_contig  & end <= end_contig,
                    "YES",
                    "NO")) %>%
    mutate(start_end_overlap = ifelse(start_overlap == "NO" & end_overlap == "NO",
                                      "NO_OVERLAP",
                                      ifelse(start_overlap == "YES" & end_overlap == "NO",
                                             "START_OVERLAP",
                                             ifelse(start_overlap == "NO" & end_overlap == "YES",
                                                    "END_OVERLAP",
                                                    ifelse(start_overlap == "YES" & end_overlap == "YES",
                                                           "BOTH_OVERLAP",
                                                           "ERROR"))))) %>% 
    select(-start_overlap,
           -end_overlap) %>%
    filter(start_end_overlap != "NO_OVERLAP",
           start_end_overlap != "ERROR",
           !is.na(start_end_overlap)) %>%
    #mutate(start = ifelse(start_end_overlap == "END_OVERLAP", start_contig, start)) %>%
    #mutate(end = ifelse(start_end_overlap == "START_OVERLAP", end_contig, end)) %>% 
    mutate(start_scaffold_new = ((((as.numeric(start) + as.numeric(start_scaffold))-1)-as.numeric(start_contig))+1)) %>%
    mutate(end_scaffols_new = ((((as.numeric(end) + as.numeric(start_scaffold))-1)-as.numeric(start_contig))+1)) %>%
    select(Scaffold,
           start = start_scaffold_new,
           end  = end_scaffols_new,
           value) ->  internal_result
  return(internal_result)
}
#convert scaffolds with fragmented contigs too original contigs positions and then convert to new scaffols coordinates
reorder_by_coord_position <- function(ARGUMENT1,ARGUMENT2,ARGUMENT3) {
ARGUMENT1 %>%
  select(Scaffold ,
         position ,
         value ) %>% 
  full_join(.,
            ARGUMENT2) %>% 
  mutate(overlap = ifelse(as.numeric(position) >= as.numeric(start_scaffold) & as.numeric(position) <= as.numeric(end_scaffold) ,"YES", "NO")) %>%
  filter(overlap == "YES") %>% 
  mutate(position  = ((((as.numeric(position) - as.numeric(start_scaffold))+1) + as.numeric(start_contig))-1 )) %>% 
  select(contig_base,
         position,
         value) %>%
  full_join(.,ARGUMENT3) %>% 
  mutate(overlap = ifelse(as.numeric(position) >= as.numeric(start_contig) & as.numeric(position) <= as.numeric(end_contig) ,"YES", "NO")) %>%
  filter(overlap == "YES") %>%
  mutate(position  = ((((as.numeric(position) + as.numeric(start_scaffold))-1) - as.numeric(start_contig))+1 )) %>%
  select(Scaffold,
         position,
         value) -> result_internal 
return(result_internal)}

plot_data <- function(DATA, #main data with columns Scaffold, position, value
                      TARGET_SCAFFOLD, #the scaffold for plot
                      Y_AXIS) { #the measurment for the plot, for example telomere signal, or LTR frequency
  return(ggplot(
        filter(DATA,
               Scaffold == TARGET_SCAFFOLD),
        aes(x=position,
            y=value)) + 
        geom_rect(
          aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold, fill = contig),
          ymin = min(filter(DATA,Scaffold == TARGET_SCAFFOLD)$value),
          ymax = max(filter(DATA,Scaffold == TARGET_SCAFFOLD)$value),
          alpha = 1,
          data = filter(coordinates,Scaffold == TARGET_SCAFFOLD)
        ) +
        geom_line() +
        #geom_point() +
        theme_classic() +
        theme(legend.position = "none",
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        ggtitle(TARGET_SCAFFOLD) +
        scale_x_continuous(name="Position (bp)", 
                           expand = c(0, 0),
                           limits = c(0,NA)) +
        scale_y_continuous(name= Y_AXIS,
                           expand = c(0, 0), 
                           limits = c(NA,NA)))
}

##########################################################Coordinations
fragmented_coordinates_modification(coordinates_file) -> coordinates
fragmented_coordinates_modification(original_coordinatess_for_EDTA) -> original_coordinates
##########################################################Linkage-map
fread(geneticmap_hap1) %>% 
  rbind(.,fread(geneticmap_hap2)) %>%
  select(Scaffold = V2,
         position = V3,
         value = V5) %>%
reorder_by_coord_position(.,original_coordinates,coordinates) -> GM

plot_data(GM,scaffold_H1,"cM")

##########################################################recombination 
fread(recomnination_hap1) %>%  
  rbind(.,fread(recomnination_hap2)) %>%
  select(Scaffold = V1,
         position = V2,
         value = V3) %>%
reorder_by_coord_position(.,original_coordinates,coordinates) %>%
  mutate(value=log10(value))-> REC

plot_data(REC,scaffold_H1,"Recombination rate")

##########################################################telomere
fread(telomere_hap1) %>% 
  rbind(fread(telomere_hap2)) %>%
  select(Scaffold = V1,
         start = V2,
         end = V3,
         value = V4) %>%
scaffold2contig_ranges(.,original_coordinates) %>%
contig2scaffold_ranges(.,coordinates) %>% 
  # gather(started,position,start:end) %>%
   mutate(value = ifelse(value == 0,0, log10(value))) %>%
  select(Scaffold,
         position = start,
         value) -> TELO

plot_data(TELO,scaffold_H1,"Telomere frequency")

##########################################################EDTA (LTR frequency)
fread(EDTA_hap1) %>% 
  rbind(.,fread(EDTA_hap2)) %>%
  filter(grepl("LTR",V2)) %>%   #(V2=="Gypsy_LTR_retrotransposon")
  select(Scaffold = V1,
         LTR_start = V3, #start of a LTR region on scaffold
         LTR_end = V4) %>%
mutate(start = floor(
  (
    (as.numeric(LTR_start)/WIND_SIZE)+(as.numeric(LTR_end)/WIND_SIZE)
    )/2)*WIND_SIZE) %>% 
  mutate(end = start + WIND_SIZE) %>%
  group_by(Scaffold,start,end) %>%
  summarize(value = n()) %>%
  ungroup() %>% 
  scaffold2contig_ranges(.,original_coordinates) %>% 
  contig2scaffold_ranges(.,coordinates) %>% 
  select(Scaffold,
         position = start,
         value) -> LTR


plot_data(LTR,scaffold_H1,"LTR frequency")






  
##########################################################Synteny
read_paf(Synteny) %>%
  as.data.frame() %>% 
  filter(mapq > 0) %>% 
  mutate(value = paste(strand,
               tname,
               tstart,
               tend,
               sep = "::")) %>%
  select(Scaffold = qname,
         start = qstart,
         end = qend,
         value)  %>%  
  scaffold2contig_ranges(.,
                         filter(original_coordinates,grepl("H2_HiC_scaffold_",Scaffold))) %>% 
    contig2scaffold_ranges(.,
                           filter(coordinates,grepl("H2_HiC_scaffold_",Scaffold)))  %>% 
  group_by(value) %>%
  mutate(newstart = min(start),
         newend = max(end)) %>% 
  select(-start,
         -end) %>%
  distinct(Scaffold,
           newstart,
           newend,
           value) %>% 
  select(Scaffold,
         start = newstart,
         end = newend,
         value) %>% 
  separate(value, 
           into = c("strand",
                    "tname",
                    "tstart",
                    "tend"),
           sep = "::",
           remove = T)  %>%
  mutate(tstart = as.numeric(tstart),
         tend = as.numeric(tend)) %>% 
  mutate(value = paste(Scaffold,
                       start,
                       end,
                       strand,
                       sep = "::")) %>% 
  select(Scaffold = tname,
         start = tstart,
         end = tend,
         value) %>% 
  scaffold2contig_ranges(.,
                         filter(original_coordinates,grepl("H1_HiC_scaffold_",Scaffold))) %>% 
  contig2scaffold_ranges(.,
                         filter(coordinates,grepl("H1_HiC_scaffold_",Scaffold))) %>% 
  group_by(value) %>%
  mutate(newstart = min(start),
         newend = max(end)) %>% 
  select(-start,
         -end) %>%
  distinct(Scaffold,
           newstart,
           newend,
           value) %>% 
  select(Scaffold,
         start = newstart,
         end = newend,
         value) %>%
  separate(value, 
           into = c("qname",
                    "qstart",
                    "qend",
                    "strand"),
           sep = "::",
           remove = T) %>%
  select(qname,
         qstart,
         qend,
         strand,
         tname = Scaffold,
         tstart = start,
         tend = end) %>%
  mutate(qstart = as.numeric(qstart),
         qend = as.numeric(qend),
         tstart = as.numeric(tstart),
         tend = as.numeric(tend)) -> alignment

alignment %>%
  filter(qname == "H2_HiC_scaffold_1",
         tname == "H1_HiC_scaffold_1")  %>% 
  mutate(tstart_new = ifelse(strand == "-",tend,tstart)) %>%
  mutate(tend_new = ifelse(strand == "-",tstart,tend)) %>% 
  select(-tstart,
         -tend) %>%
  rename(tstart = tstart_new,
         tend = tend_new) %>%
ggplot(.) + 
  geom_rect(
    aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold,fill=contig),
    #fill = NA,
    #colour = "black",
    ymin = 0,
    ymax = max(max(alignment$qend),max(alignment$qstart)),
    alpha = 1,
    data = filter(coordinates,Scaffold == scaffold_H1)) +
  geom_segment(aes(x = qstart, 
                   xend = qend, 
                   y = tstart, 
                   yend = tend
  )) +
  geom_rect(
    aes(NULL,NULL,ymin = start_scaffold, ymax = end_scaffold),
    fill = NA,
    colour = "black",
    xmin = 0,
    xmax = max(max(ali1$tend),max(ali1$tstart)),,
    alpha = 0.2,
    data = filter(coor,Scaffold == scaffold_H2)) +
  

###
  ggplot(
    filter(DATA,
           Scaffold == TARGET_SCAFFOLD),
    aes(x=position,
        y=value)) + 
  geom_rect(
    aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold, fill = contig),
    ymin = min(filter(DATA,Scaffold == TARGET_SCAFFOLD)$value),
    ymax = max(filter(DATA,Scaffold == TARGET_SCAFFOLD)$value),
    alpha = 1,
    data = filter(coordinates,Scaffold == TARGET_SCAFFOLD)
  ) +
  geom_line() +
  #geom_point() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle(TARGET_SCAFFOLD) +
  scale_x_continuous(name="Position (bp)", 
                     expand = c(0, 0),
                     limits = c(0,NA)) +
  scale_y_continuous(name= Y_AXIS,
                     expand = c(0, 0), 
                     limits = c(NA,NA))) 
###


read_paf(Synteny) %>%
  as.data.frame() %>% 
  filter(mapq > 0) %>%
  mutate(tstart_new = ifelse(strand == "-",tend,tstart)) %>%
  mutate(tend_new = ifelse(strand == "-",tstart,tend)) %>% 
  select(-tstart,
         -tend) %>%
  rename(tstart = tstart_new,
         tend = tend_new) %>%
  filter(qname=="H2_HiC_scaffold_1",
         tname=="H1_HiC_scaffold_1") -> b %>% 
ggplot(.) + 
  geom_segment(aes(x = qstart, 
                   xend = qend, 
                   y = tstart, 
                   yend = tend
                   #color=factor(tname)
  ))
  # geom_rect(
  #   aes(NULL,NULL,xmin = tstart, xmax = tend),
  #   #fill = NA,
  #   #colour = "black",
  #   ymin = 0,
  #   ymax = max(max(alignment$qend),max(alignment$qstart)),
  #   alpha = 1,
  #   data = filter(coor,Scaffold == scaffold_H1)) +
  # geom_rect(
  #   aes(NULL,NULL,ymin = start, ymax = end),
  #   fill = NA,
  #   colour = "black",
  #   xmin = 0,
  #   xmax = max(max(ali1$tend),max(ali1$tstart)),,
  #   alpha = 0.2,
  #   data = filter(coor,Scaffold == scaffold_H2)) +
  geom_segment(aes(x = tstart, 
                   xend = tend, 
                   y = qstart, 
                   yend = qend
                   #color=factor(tname)
  )) + 
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = scaffold_H1, 
       y = scaffold_H2) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

    
    
  left_join(.,
            select(original_coordinates,
                   scaffold = Scaffold,qname, start_scaff=start))  %>% 
  original_coordinates %>% glimpse
  
  ##########################################################DEPTH
  fread(DEPTH_hap1) %>% 
    rbind(fread(DEPTH_hap2))  %>% 
    select(Scaffold = V1,
           start = V2,
           end = V3,
           value = V4) %>% 
    scaffold2contig_ranges(.,original_coordinates) %>%
    contig2scaffold_ranges(.,coordinates) -> DEPTH_hap1
  
  fread("/Users/mojtabajahani/Downloads/AGA10_r0/new_30_june_2022/depth/depth.bed") -> DEPTH 
  mutate(depth=ifelse(depth == 0,
                      0,
                      (ceiling(depth/5))*5))
  
  fread(DEPTH_hap1) %>% 
    rbind(fread(DEPTH_hap2))  %>% 
    select(Scaffold = V1,
           start = V2,
           end = V3,
           value = V4) %>%
    filter(Scaffold=="H1_HiC_scaffold_1") -> a
  mean(a$value)-> HOMO
  HOMO/2 -> HETRO
  HETRO/2 -> RANGE
  
  a %>%
    mutate(new_value=ifelse(value < RANGE,
                            RANGE/2,
                            ifelse(between(value,RANGE,RANGE+HETRO),
                                   HETRO,
                                   ifelse(between(value,RANGE+HETRO,RANGE+HOMO),
                                          HOMO,
                                          HOMO+(2*RANGE))))) %>% ggplot(.,aes(start,new_value)) +  geom_point() -> ab
  ggplotly(ab)
  
  mutate(value=ifelse(value == 0,
                      0,
                      (ceiling(value/5))*5)) %>% 
    ggplot(.,aes(start,value)) +  geom_line() -> a
  ggplotly(a) 
  ggplot(
    filter(DEPTH,
           Scaffold == "H1_HiC_scaffold_1"),
    aes(x = start,
        y = value)) + 
    geom_rect(
      aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold, fill = contig),
      ymin = min(filter(LTR,Scaffold == scaffold_H1)$value),
      ymax = max(filter(LTR,Scaffold == scaffold_H1)$value),
      alpha = 1,
      data = filter(coordinates,Scaffold == scaffold_H1)
    ) +
    geom_line() +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(scaffold_H1)  +
    scale_x_continuous(name="Position (bp)", 
                       expand = c(0, 0),
                       limits = c(0,NA)) +
    scale_y_continuous(name="Telomere frequency",
                       expand = c(0, 0), 
                       limits = c(NA,NA)
  
########################
fread(telomere_hap1) %>% 
  rbind(fread(telomere_hap2)) %>%
  select(Scaffold = V1,
         start = V2,
         end = V3,
         value = V4)-> b
a -> TELO


ggplot(
  filter(LTR,
         Scaffold == scaffold_H1),
  aes(x=start,
      y=value)) + 
  geom_rect(
    aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold, fill = contig),
    ymin = 0,
    ymax = max(filter(LTR,Scaffold == scaffold_H1)$value),
    alpha = 1,
    data = filter(coordinates,Scaffold == scaffold_H1)
  ) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle(scaffold_H1) +
  scale_x_continuous(name="Position (bp)", 
                     expand = c(0, 0),
                     limits = c(0,NA)) +
  scale_y_continuous(name="Telomere frequency",
                     expand = c(0, 0), 
                     limits = c(NA,NA)) 

set.seed(1234)
df <- data.frame(
  sex=factor(rep(c("F", "M"), each=200)),
  weight=round(c(rnorm(200, mean=55, sd=5),
                 rnorm(200, mean=65, sd=5)))
)
head(df)
dat <- with(density(df$weight), data.frame(x, y))
dat %>% head

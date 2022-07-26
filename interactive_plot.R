library(plotly)
library(data.table)
library(tidyverse)
library(pafr)
library(htmlwidgets)
options(scipen = 9999)
WIND_SIZE=200000 #EDTA window size
args = commandArgs(trailingOnly = TRUE)
##########################################################
original_coordinates_file  <- args[1]
coordinates_file <- args[2]

geneticmap_hap1 <- args[3]
geneticmap_hap2 <- args[4]

recomnination_hap1 <- args[5]
recomnination_hap2 <- args[6]

telomere_hap1 <- args[7]
telomere_hap2 <- args[8]

EDTA_hap1 <- args[9]
EDTA_hap2 <- args[10]

Synteny <- args[11]
SAVE_DIR <- args[12]
##########################################################Functions
#calculate the start and end position of fragmented contigs and their coordinates in scaffolds
fragmented_coordinates_modification <- function(ARGUMENT1) {
  fread(ARGUMENT1) %>% 
    select(Scaffold = V1,
           contig = V4,
           start = V2,
           end = V3,
           orientation = V5) %>% 
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
             end_contig,
             orientation) -> FRAGMENTED_SE
    
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
             end_contig,
             orientation) %>%
      rbind(.,FRAGMENTED_SE) -> coordinates_internal
    
    
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
             end_contig,
             orientation) -> coordinates_internal
  }
  
  return(coordinates_internal)
  
  rm(FRAGMENTS,
     FRAGMENTED_SE,
     contig_coor)
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
           value,
           orientation_original = orientation) -> internal_result
  return(internal_result)
}

#original contigs to new fragmented contig scaffold position

contig2scaffold_ranges <- function(ARGUMENT1,#result of scaffold2contig
                                   ARGUMENT2) { #new coordinates
  ARGUMENT1 %>%
    select(contig_base,
           start,
           end,
           value,
           orientation_original)  %>% 
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
    mutate(sign = as.character(orientation_original == orientation)) %>% 
    #mutate(start = ifelse(start_end_overlap == "END_OVERLAP", start_contig, start)) %>%
    #mutate(end = ifelse(start_end_overlap == "START_OVERLAP", end_contig, end)) %>% 
    mutate(start_scaffold_new = ifelse(sign == "TRUE",
                                       (((as.numeric(start) - as.numeric(start_contig))-1)+ as.numeric(start_scaffold)+1),
                                       ((as.numeric(end_contig) - as.numeric(start))-1) + as.numeric(start_scaffold)+1)
    ) %>% 
    mutate(end_scaffols_new = ifelse(sign == "TRUE",
                                     (((as.numeric(end) - as.numeric(start_contig))-1)+ as.numeric(start_scaffold)+1),
                                     ((as.numeric(end_contig) - as.numeric(end))-1) + as.numeric(start_scaffold)+1)
    ) %>% 
    select(Scaffold,
           start = start_scaffold_new,
           end  = end_scaffols_new,
           value,
           sign) -> internal_result
  return(internal_result)}

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
           value,
           orientation_original = orientation) %>% 
    full_join(.,ARGUMENT3) %>% 
    mutate(overlap = ifelse(as.numeric(position) >= as.numeric(start_contig) & as.numeric(position) <= as.numeric(end_contig) ,"YES", "NO")) %>%
    filter(overlap == "YES")  %>% 
    mutate(position  = ifelse(orientation_original == orientation ,
                              (((as.numeric(position) - as.numeric(start_contig))-1)+ as.numeric(start_scaffold)+1),
                              ((as.numeric(end_contig) - as.numeric(position))-1) + as.numeric(start_scaffold)+1)
    ) %>%
    select(Scaffold,
           position,
           value)  -> result_internal 
  return(result_internal)}

plot_data <- function(DATA, #main data with columns Scaffold, position, value
                      TARGET_SCAFFOLD, #the scaffold for plot
                      Y_AXIS,#the measurment for the plot, for example telomere signal, or LTR frequency
                      POINT = "NO") { #geom_pont
  if (POINT=="YES") {
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
        geom_point() +
        theme_classic() +
        theme(legend.position = "none"
        ) +
        scale_x_continuous(name="Position (bp)", 
                           expand = c(0, 0),
                           limits = c(0,NA)) +
        scale_y_continuous(name= Y_AXIS,
                           expand = c(0, 0), 
                           limits = c(NA,NA)))
  } else { 
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
        theme_classic() +
        theme(legend.position = "none"
        ) +
        scale_x_continuous(name="Position (bp)", 
                           expand = c(0, 0),
                           limits = c(0,NA)) +
        scale_y_continuous(name= Y_AXIS,
                           expand = c(0, 0), 
                           limits = c(NA,NA)))}
  
}

##########################################################Coordinates
fragmented_coordinates_modification(coordinates_file) -> coordinates #calculate the start and end position of fragmented contigs and their coordinates in scaffolds
fragmented_coordinates_modification(original_coordinates_file) -> original_coordinates
##########################################################Linkage-map
fread(geneticmap_hap1) %>% 
  rbind(.,fread(geneticmap_hap2)) %>%
  select(Scaffold = V2,
         position = V3,
         value = V5) %>% 
  reorder_by_coord_position(.,original_coordinates,coordinates) -> GM
##########################################################recombination 
fread(recomnination_hap1) %>%  
  rbind(.,fread(recomnination_hap2)) %>%
  select(Scaffold = V1,
         position = V2,
         value = V3) %>%
  reorder_by_coord_position(.,original_coordinates,coordinates) %>%
  mutate(value=log10(value)) -> REC
##########################################################telomere
fread(telomere_hap1) %>% 
  rbind(fread(telomere_hap2)) %>%
  select(Scaffold = V1,
         start = V2,
         end = V3,
         value = V4) %>% 
  scaffold2contig_ranges(.,original_coordinates) %>% 
  contig2scaffold_ranges(.,coordinates) %>% 
  mutate(value = ifelse(value == 0,0, log10(value))) %>%
  select(Scaffold,
         position = start,
         value) -> TELO
##########################################################EDTA (LTR frequency)
fread(EDTA_hap1) %>% 
  rbind(.,fread(EDTA_hap2)) %>%
  filter(grepl("LTR",V2)) %>%  #(V2=="Gypsy_LTR_retrotransposon")
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
##########################################################Synteny
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
         tend = new_tend)-> SYNTENY
###########################################################################################################
################################################DRAWING PLOTS##############################################
###########################################################################################################
for (CHR in seq(10)) {
  paste0("H1_HiC_scaffold_",CHR) -> scaffold_H1
  paste0("H2_HiC_scaffold_",CHR) -> scaffold_H2
  ggplotly(plot_data(GM,scaffold_H1,"cM","YES"),dynamicTicks=T) %>% layout(annotations = list(x = 0.5 , y = 1.2,text = scaffold_H1 ,showarrow = F, xref='paper', yref='paper')) -> Geneticmap_plot_hap1 #layout(title = scaffold_H1)
  ggplotly(plot_data(GM,scaffold_H2,"","YES"),dynamicTicks=T) %>% layout(annotations = list(x = 0.5 , y = 1.2,text = scaffold_H2, showarrow = F, xref='paper', yref='paper')) -> Geneticmap_plot_hap2
  ggplotly(plot_data(REC,scaffold_H1,"Recombination"),dynamicTicks=T) -> Recombinationrate_plot_hap1 
  ggplotly(plot_data(REC,scaffold_H2,""),dynamicTicks=T) -> Recombinationrate_plot_hap2
  ggplotly(plot_data(TELO,scaffold_H1,"Telomere"),dynamicTicks=T) -> Telomerefrequency_plot_hap1
  ggplotly(plot_data(TELO,scaffold_H2,""),dynamicTicks=T) -> Telomerefrequency_plot_hap2
  ggplotly(plot_data(LTR,scaffold_H1,"LTR"),dynamicTicks=T) -> LTRfrequency_plot_hap1
  ggplotly(plot_data(LTR,scaffold_H2,""),dynamicTicks=T) -> LTRfrequency_plot_hap2
  
  
  SYNTENY %>%
    filter(qname == scaffold_H2,
           tname == scaffold_H1) -> alignment
    
  alignment %>%
    ggplot(.) + 
    geom_rect(
      aes(NULL,NULL,ymin = start_scaffold, ymax = end_scaffold), 
      fill = NA,
      colour = "black",
      xmin = 0,
      xmax = max(max(alignment$tend),max(alignment$tstart)),
      alpha = 1,
      data = filter(coordinates,Scaffold == scaffold_H2)) +
    geom_rect(
      aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold,fill=contig),
      ymin = 0,
      ymax = max(max(alignment$qend),max(alignment$qstart)), 
      alpha = 0.95,
      data = filter(coordinates,Scaffold == scaffold_H1)) + 
    theme(legend.position = "none") +
    geom_segment(aes(x = tstart, 
                     xend = tend, 
                     y = qstart, 
                     yend = qend
    )) +
    scale_x_continuous(name ="Position (bp)", 
                       expand = c(0, 0),
                       limits = c(0,NA)) +
    scale_y_continuous(name = "Position",
                       expand = c(0, 0), 
                       limits = c(NA,NA)) -> Synteny_plot_hap1
  ggplotly(Synteny_plot_hap1,dynamicTicks=T) -> Synteny_plot_hap1
  
  alignment %>%
    filter(qname == scaffold_H2,
           tname == scaffold_H1) %>%
    ggplot(.) + 
    geom_rect(
      aes(NULL,NULL,ymin = start_scaffold, ymax = end_scaffold),
      fill = NA,
      colour = "black",
      xmin = 0,
      xmax = max(max(alignment$qend),max(alignment$qstart)),
      alpha = 1,
      data = filter(coordinates,Scaffold == scaffold_H1)) +
    geom_rect(
      aes(NULL,NULL,xmin = start_scaffold, xmax = end_scaffold,fill=contig),
      ymin = 0,
      ymax = max(max(alignment$tend),max(alignment$tstart)),
      alpha = 0.95,
      data = filter(coordinates,Scaffold == scaffold_H2)) +
    theme(legend.position = "none") +
    geom_segment(aes(x = qstart, 
                     xend = qend, 
                     y = tstart, 
                     yend = tend
    )) +
    scale_x_continuous(name ="Position", 
                       expand = c(0, 0),
                       limits = c(0,NA)) +
    scale_y_continuous(name = "",
                       expand = c(0, 0), 
                       limits = c(NA,NA)) -> Synteny_plot_hap2
  ggplotly(Synteny_plot_hap2,dynamicTicks=T) -> Synteny_plot_hap2
  
  
  subplot(Geneticmap_plot_hap1,
          Geneticmap_plot_hap2,
          Recombinationrate_plot_hap1,
          Recombinationrate_plot_hap2,
          Telomerefrequency_plot_hap1,
          Telomerefrequency_plot_hap2,
          LTRfrequency_plot_hap1,
          LTRfrequency_plot_hap2,
          Synteny_plot_hap1,
          Synteny_plot_hap2,
          nrows = 5, 
          margin = 0.01, 
          heights = c(0.15, 0.15,0.15, 0.15,0.4),
          shareX = T,
          titleY = T) -> FINAL_PLOT
  
  saveWidget(FINAL_PLOT, paste0(SAVE_DIR,"/Chr",CHR,".html"), selfcontained = F, libdir = "lib")
}

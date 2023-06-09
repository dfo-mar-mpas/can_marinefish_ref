###Summary plot of the database -- using outputs from taxonomic process.R 

#load libraries----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(sf)
library(rnaturalearth)

#load data from taxonomic_process.R
load("output/ocean_intersection.RData")
load("output/ocean_intersection_lg.RData")
load("output/sp_fishbase.RData")
load("output/taxonomy_wide.RData")

#add in the subspecies
sp_fishbase 

#compile a dataframe 
ocean_df1 <- do.call("rbind",ocean_intersection_lg)%>%
  rbind(.,do.call("rbind",ocean_intersection))%>%
  arrange(species,ocean)%>%
  left_join(.,sp_fishbase)

#fix subspecies
ocean_df2 <- rbind(ocean_df1,
                   ocean_df1[ocean_df1$species == "Esox americanus",], #two subspecies for Exox americanus
                   ocean_df1[ocean_df1$species == "Oncorhynchus clarkii",], #three subspecies for Onorhynchus clarkii
                   ocean_df1[ocean_df1$species == "Oncorhynchus clarkii",])

ocean_df2$subspecies <- NA

ocean_df2[ocean_df2$species=="Esox americanus","subspecies"] <- c("Esox americanus americanus","Esox americanus vermiculatus")
ocean_df2[ocean_df2$species=="Oncorhynchus clarkii","subspecies"] <- c("Oncorhynchus clarkii bouvieri","Oncorhynchus clarkii clarkii","Oncorhynchus clarkii lewisi")

#clean it up
out_df <- out_df%>%
  mutate(spp = case_when(is.na(subspecies)~species,
                         !is.na(subspecies)~subspecies))

#the process did pick up some duplicates with on 'itis' and one 'gbif' registry. 
dup_sp <- out_df%>%
  pull(spp)%>%
  table()%>%
  data.frame()%>%
  arrange(-Freq)%>%
  rename(spp=1)%>%
  filter(Freq>1)%>%
  pull(spp)%>%
  as.character()

tax_df <- out_df%>%
  filter(!spp %in% dup_sp | (spp %in% dup_sp & db=="itis")) #keep the itis data for the duplicated species.

ocean_df <- ocean_df2%>%
  mutate(spp = case_when(is.na(subspecies)~species,
                         !is.na(subspecies)~subspecies))%>%
  dplyr::select(spp,ocean,dist)%>%
  left_join(.,tax_df,by="spp")

ocean_table <- <- 
  
  
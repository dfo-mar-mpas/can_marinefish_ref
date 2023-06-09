###Summary plot of the database -- using outputs from taxonomic process.R 

#load libraries----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(sf)
library(rnaturalearth)
library(eulerr)

#load data from taxonomic_process.R
load("output/ocean_intersection.RData")
load("output/ocean_intersection_lg.RData")
load("output/sp_fishbase.RData")
load("output/taxonomy_wide.RData")

#plot order
ocean_ord <- c("Atlantic","Arctic","Pacific")

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

#add in the subspecies
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

ocean_table <- table(ocean_df$spp,ocean_df$ocean)%>%
               data.frame()%>%
               mutate(Freq=as.logical(Freq))%>%
               spread(key=Var2,value=Freq)%>%
               dplyr::select(-Var1)

png(filename = "inst/plot_venn.png",height=6,width=6,units="in",res=300)               
  plot(euler(ocean_table, shape = "ellipse"), quantities = TRUE,fills = c("lightskyblue", "royalblue4", "springgreen4"))
dev.off()

ocean_table2 <- ocean_df%>%
                mutate(ocean=factor(ocean,levels=ocean_ord))%>%
                arrange(ocean)%>%
                group_by(spp)%>%
                summarise(ocean=paste(ocean,collapse="-"))%>%
                ungroup()%>%
                data.frame()%>%
                pull(ocean)%>%
                table()%>%
                data.frame()%>%
                rename(ocean=1,count=2)%>%
                arrange(-count)

ocean_table3 <- ocean_df%>%
                group_by(ocean)%>%
                summarise(count=n())%>%
                ungroup()%>%
                data.frame()%>%
                mutate(ocean=factor(ocean,levels=ocean_ord))%>%
                arrange(ocean)

draw.triple.venn(area1=ocean_table3%>%filter(ocean=="Atlantic")%>%pull(count),
                 area2=ocean_table3%>%filter(ocean=="Arctic")%>%pull(count),
                 area3=ocean_table3%>%filter(ocean=="Pacific")%>%pull(count),
                 n12=ocean_table2%>%filter(ocean=="Atlantic-Arctic")%>%pull(count)-ocean_table2%>%filter(ocean=="Atlantic-Arctic-Pacific")%>%pull(count),
                 n23=ocean_table2%>%filter(ocean=="Arctic-Pacific")%>%pull(count)-ocean_table2%>%filter(ocean=="Atlantic-Arctic-Pacific")%>%pull(count),
                 n13=ocean_table2%>%filter(ocean=="Atlantic-Pacific")%>%pull(count)-ocean_table2%>%filter(ocean=="Atlantic-Arctic-Pacific")%>%pull(count),
                 n123=ocean_table2%>%filter(ocean=="Atlantic-Arctic-Pacific")%>%pull(count))

n123=ocean_table2%>%filter(ocean=="Atlantic-Arctic-Pacific")%>%pull(count)
  
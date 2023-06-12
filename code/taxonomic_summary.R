###Summary plot of the database -- using outputs from taxonomic process.R 

#load libraries----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(sf)
library(rnaturalearth)
library(eulerr)
library(gridExtra)

sf_use_s2(FALSE)

#load data from taxonomic_process.R
load("output/ocean_df.RData")
load("output/taxonomy_wide.RData")

#projections -----------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0" #works for the Scotian Shelf
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#plot order --------
ocean_ord <- c("Atlantic","Arctic","Pacific")

#polygons ----------
bioregions <- read_sf("data/shapefiles/DFO_Marine_Bioregions_Clipped_1M_CAEAC_2012_05_31.shp")%>%
              st_transform(CanProj)

oceans <- bioregions%>%
          filter(REGION %in% ocean_ord)%>%
          group_by(REGION)%>%
          summarise(geometry=st_union(geometry))%>%
          ungroup()%>%
          st_as_sf()%>%
          st_transform(CanProj)%>%
          mutate(region=factor(REGION,levels=ocean_ord))%>% #sort
          arrange(region)%>%
          dplyr::select(-region)

can_bounding <- bioregions%>%st_bbox()

can_bounding_polygon <- can_bounding%>%st_as_sfc()%>%st_as_sf() #convert to an sf polygon

ggplot()+geom_sf(data=oceans)+geom_sf(data=can_bounding_polygon,fill=NA)

#create a basemap for plotting MPAs
basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(latlong)%>%
                   st_as_sf()%>%
                   mutate(country="Canada"),
                 ne_states(country = "United States of America",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(latlong)%>%
                   st_as_sf()%>%
                   mutate(country="USA"),
                 ne_states(country = "Greenland",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(latlong)%>%
                   st_as_sf()%>%
                   mutate(country="USA"))%>%
  st_union()%>%
  st_as_sf()%>%
  st_transform(CanProj)%>%
  st_intersection(.,can_bounding_polygon)

Canada <-  ne_states(country = "Canada",returnclass = "sf")%>%
            dplyr::select(name_en,geometry)%>%
            st_as_sf()%>%
            st_union()%>%
            st_transform(CanProj)%>%
            st_as_sf()%>%
            mutate(country="Canada")

#Data for venn diagram -----

ocean_table1 <- ocean_df%>%
               mutate(Atlantic = ifelse(grepl("Atlantic",ocean),TRUE,FALSE),
                      Arctic = ifelse(grepl("Arctic",ocean),TRUE,FALSE),
                      Pacific = ifelse(grepl("Pacific",ocean),TRUE,FALSE),
                      Freshwater = ifelse(env == "Freshwater",TRUE,FALSE))%>%
               dplyr::select(species,Atlantic,Arctic,Pacific,Freshwater)%>%
               gather(key='Ocean',value='logic',Atlantic:Freshwater)%>%
               filter(logic)%>%
               data.frame()

ocean_table <- table(ocean_table1$species,ocean_table1$Ocean)%>%
               data.frame()%>%
                mutate(Freq=as.logical(Freq))%>%
                spread(key=Var2,value=Freq)%>%
                dplyr::select(-Var1)


#make the venn plot
set.seed(23) #so you get the same bubble orientation each time. 

venn.plot <- plot(euler(ocean_table, shape = "ellipse"), quantities = TRUE,fills = c("lightskyblue", "royalblue4", "forestgreen","seagreen2"))

png(filename = "inst/plot_venn.png",height=6,width=6,units="in",res=300) 
venn.plot 
dev.off()

## Companion map

p1 <- ggplot()+
  geom_sf(data=basemap,fill="white")+
  geom_sf(data=Canada,fill="forestgreen")+
  geom_sf(data=oceans,aes(fill=REGION))+
  scale_fill_manual(values=c("lightskyblue", "royalblue4", "seagreen2"))+
  theme_bw()+
  coord_sf(expand=0,xlim=can_bounding[c(1,3)],ylim=can_bounding[c(2,4)])+
  theme(legend.position="none")

ggsave("inst/canada_map.png",p1,height=5,width=6,units="in",dpi=300)

#combination for github readme
png(filename = "inst/combination.png",height=12,width=6,units="in",res=300) 
 grid.arrange(test,p1,nrow=2)
dev.off()
  
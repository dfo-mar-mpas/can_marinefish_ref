## extended search of the Canadian Adjacent waters (US North Atlantic, Pacific, and Alaska)

library(dplyr)
library(sf)
library(robis)
library(mregions)
library(broom)
library(rnaturalearth)
library(worrms)
library(tidyr)
library(patchwork)

sf_use_s2(FALSE)

#load the classification 
source("code/worms_classify.R")

#projections -----------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the Canadian species list
load("output/can_spec_list.RData")

#load the oceans polygon generated in taxonomic_process.R
load("data/oceans.RData")
oceans <- oceans%>%
          st_transform(latlong)%>%
          mutate(areaid = case_when(REGION == "Atlantic" ~ 272, #Matching the robis mregions areaids
                                    REGION == "Arctic" ~ 265,
                                    REGION == "Pacific" ~ 274))

#basemaps 
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
  st_transform(CanProj)

Canada <-  ne_states(country = "Canada",returnclass = "sf")%>%
  dplyr::select(name_en,geometry)%>%
  st_as_sf()%>%
  st_union()%>%
  st_transform(CanProj)%>%
  st_as_sf()%>%
  mutate(country="Canada")

#Search for species in Canadian adjacent Atlantic and pacific waters
us_alaska <- checklist(areaid = 265)%>%as_tibble()%>%mutate(ocean="Alaska",areaid = 265)
us_north_pacific <- checklist(areaid = 274)%>%as_tibble()%>%mutate(ocean="North Pacific",areaid=274)
us_north_atlantic <- checklist(areaid = 272)%>%as_tibble()%>%mutate(ocean="North Atlantic",areaid=272)

us_names <- intersect(names(us_alaska),names(us_north_pacific))%>%intersect(.,names(us_north_atlantic)) #they have some differences

us_taxa <- rbind(us_alaska%>%dplyr::select(all_of(us_names)),
                 us_north_pacific%>%dplyr::select(all_of(us_names)),
                 us_north_atlantic%>%dplyr::select(all_of(us_names)))

#save interim output
save(us_taxa,file="output/obis_checklist_US.RData")

#format data 
us_df <- us_taxa%>%
  filter(taxonomicStatus == 'accepted',#keep if it has a valid iD
         as.logical(is_marine)|as.logical(is_brackish), #keep if it is marine or brackish
         superclass=="Actinopteri" | subclass=="Teleostei" | class=="Teleostei",
         !is.na(species))%>% #keep the bony fishes
  dplyr::select(scientificName,areaid,taxonID,bold_id,ncbi_id,kingdom,phylum,class,order,family,genus,species)%>%
  rename(AphiaID = taxonID)

#identify the new taxa not within the databases
us_missing <- setdiff(us_df%>%pull(AphiaID)%>%unique(),can_spec_list%>%pull(AphiaID))

#create some buffering polygons to speed up the process
prox_thresh <- 1500 # 500 km of the coast. 

pacific_buffer <- oceans%>%
                  filter(REGION=="Pacific")%>%
                  mutate(areaid=274)%>%
                  st_transform(CanProj)%>%
                  st_buffer(prox_thresh*1000)%>%
                  st_simplify()

atlantic_buffer <- oceans%>%
                    filter(REGION=="Atlantic")%>%
                    mutate(areaid=272)%>%
                    st_transform(CanProj)%>%
                    st_buffer(prox_thresh*1000)%>%
                    st_simplify()

arctic_buffer <- oceans%>%
                filter(REGION=="Arctic")%>%
                mutate(areaid=265)%>%
                st_transform(CanProj)%>%
                st_buffer(prox_thresh*1000)%>%
                st_difference(.,atlantic_buffer)#since migrants will come from the Canadian region and the North Atlantic, the search doesn't need to include this area. 
                st_simplify()

buffers <- rbind(pacific_buffer,atlantic_buffer,arctic_buffer)%>%st_transform(latlong)

#plot it to visualize the buffers

pac_plot <- ggplot()+
            geom_sf(data=pacific_buffer,fill="seagreen2",alpha=0.5)+
            geom_sf(data=basemap)+
            geom_sf(data=oceans%>%filter(REGION=="Pacific")%>%st_transform(CanProj),fill="seagreen2")+
            coord_sf(xlim = st_bbox(pacific_buffer)[c(1,3)],ylim= st_bbox(pacific_buffer)[c(2,4)],expand=0)+
            theme_bw()+
            labs(title=paste0("Pacific buffer ",prox_thresh,"km"))

atl_plot <- ggplot()+
            geom_sf(data=atlantic_buffer,fill="royalblue4",alpha=0.5)+
            geom_sf(data=basemap)+
            geom_sf(data=oceans%>%filter(REGION=="Atlantic")%>%st_transform(CanProj),fill="royalblue4")+
            coord_sf(xlim = st_bbox(atlantic_buffer)[c(1,3)],ylim= st_bbox(atlantic_buffer)[c(2,4)],expand=0)+
            theme_bw()+
            labs(title=paste0("Atlantic buffer ",prox_thresh,"km"))

arc_plot <- ggplot()+
            geom_sf(data=arctic_buffer,fill="lightskyblue",alpha=0.5)+
            geom_sf(data=basemap)+
            geom_sf(data=oceans%>%filter(REGION=="Arctic")%>%st_transform(CanProj),fill="lightskyblue")+
            coord_sf(xlim = st_bbox(arctic_buffer)[c(1,3)],ylim= st_bbox(arctic_buffer)[c(2,4)],expand=0)+
            theme_bw()+
            labs(title=paste0("Arctic buffer ",prox_thresh,"km"))

buffer_plot <- pac_plot + arc_plot + atl_plot + plot_layout(ncol=3)

ggsave("output/buffer_plot.png",buffer_plot,width=12,height=4,units="in",dpi=300)
            
#run obis searches - and extract staight line distances based on a proximity threshold

areas <- data.frame(region=c("Alaska","North Pacific","North Atlantic"),
                    areaid=c(265,274,272),
                    ocean=c("Arctic","Pacific","Atlantic"))#corresponding to the search buffers. 

us_sp <- us_df%>%
  filter(AphiaID %in% us_missing)%>%
  dplyr::select(scientificName,areaid)

sp_list <- us_sp%>%pull(scientificName)%>%unique()

us_extracts <- list()

for(i in sp_list){
  
  temp <- us_sp%>%filter(scientificName==i)
  
  for(j in unique(temp$areaid)){
  
   message(paste0("Working on ",i," ",which(i == sp_list),"/",length(sp_list)," species, in ",areas%>%filter(areaid==j)%>%pull(region)," - ",which(j == unique(temp$areaid)),"/",length(unique(temp$areaid))," unique regions."))
  
  #check to see if the file has been extracted  
  if(!file.exists(paste0("output/robis_records/us_records/robis","_AreaID-",j,"_",i,".RData"))){
    
   temp2 <- robis::occurrence(taxonid=us_df%>%
                                      filter(scientificName==i)%>%
                                      distinct(AphiaID)%>%
                                      pull(AphiaID),
                              areaid = j)
   
   #save the obis extraction
   save(temp2,file=paste0("output/robis_records/us_records/robis","_AreaID-",j,"_",i,".RData"))
   
    } else {load(paste0("output/robis_records/us_records/robis","_AreaID-",j,"_",i,".RData"))} #load it if it already exists
   
    message("intersection ....")
    
    if(!is.null(temp2)){
   
   temp3 <- temp2%>%
            dplyr::select(aphiaID,species,decimalLongitude,decimalLatitude)%>%
            st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
            st_intersection(.,buffers%>%filter(areaid == j))%>%
            suppressWarnings()%>%
            suppressMessages()
   
       #now find the distances 
       if(nrow(temp3)>0){
         
         message("calculating distance to eez ....")
         
         temp4 <- temp3%>%
                  mutate(distance=as.numeric(st_distance(.,oceans%>%filter(areaid==j)))/1000)%>%
                  distinct(distance,.keep_all = TRUE)%>%
                  filter(distance==min(distance,na.rm=T))
         
         out <- temp4%>%
                data.frame()%>%
                dplyr::select(aphiaID,species,REGION,areaid,distance)
         
       } 
   
   
        if(nrow(temp3) == 0){out <- data.frame(aphiaID = us_df%>%filter(scientificName==i)%>%distinct(AphiaID)%>%pull(AphiaID),
                                species=i,
                                REGION=buffers%>%filter(areaid == j)%>%pull(REGION),
                                areaid=j,
                                distance=NA)
                            }
    
    }
    
    if(is.null(temp2)){out <- data.frame(aphiaID = us_df%>%filter(scientificName==i)%>%distinct(AphiaID)%>%pull(AphiaID),
                            species=i,
                            REGION=buffers%>%filter(areaid == j)%>%pull(REGION),
                            areaid=j,
                            distance=NA)}
            
  us_extracts[[paste(i,j,sep="-")]] <- out   
  
  rm(out)
  
  }#j loop end
  #interim save incase it crashes
  save(us_extracts,file="output/us_extracts.RData") #will overwrite each loop
}# i loop end

## do the worrms extraction for the environment

worms_classification_us <- list()

for(i in 1:length(sp_list)){
  
  message(paste0("Working on ",sp_list[i]," ",i,"/",length(sp_list)))
  
  worms_classification_us[[sp_list[i]]] <- worms_classify(sp_list[i])
  
}

save(worms_classification_us,file="output/worms_classification_us.RData")


#compile dataframe

us_extracts_df <- NULL

for(i in ocean_ord){

us_extracts_df <- us_extracts_df%>%
                  rbind(.,
                    do.call("rbind",us_extracts)%>%
                    filter(REGION == i)%>%
                    left_join(us_df%>%
                                left_join(.,areas%>%rename(REGION=ocean)%>%dplyr::select(areaid,REGION))%>%
                                filter(REGION==i)%>%
                                rename(aphiaID=AphiaID)%>%
                                dplyr::select(aphiaID,bold_id,ncbi_id,kingdom,phylum,class,order,family,genus))%>%
                    dplyr::select(REGION,areaid,distance,aphiaID,bold_id,ncbi_id,kingdom,phylum,class,order,family,genus,species))
}

row.names(us_extracts_df) <- NULL

us_extracts_df <- us_extracts_df%>%filter(!is.na(distance)) #na values didn't fall within the buffered areas. 

#save formatted output
save(us_extracts_df,file="output/us_extracts_df.RData")

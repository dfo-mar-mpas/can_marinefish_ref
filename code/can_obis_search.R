## code to evaluate all marine fish species identified within Canada based on the 
## original database (He et al. 2022/2023) species identified within the Canadian Register of Marine Species (CaRMS)
## and OBIS

##load libraries ------
library(dplyr)
library(sf)
library(robis)
library(mregions)
library(broom)
library(rnaturalearth)
library(worrms)
library(tidyr)

sf_use_s2(FALSE)

source('code/worms_classify.R')

#projections -----------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the data extracted from the original database and validated using worms -----
load("output/taxonomy_wide.RData")

#load polygons -------------
    bioregions <- read_sf("data/shapefiles/DFO_Marine_Bioregions_Clipped_1M_CAEAC_2012_05_31.shp")
    
    ocean_ord <- c("Atlantic","Arctic","Pacific")

  #search extent polygons 
      bb_atlantic <- bioregions%>%
        filter(REGION == "Atlantic")%>%
        st_bbox()%>%
        st_as_sfc()%>%
        st_transform(CanProj)%>% #area based projection (units m)
        st_buffer(1000 * 1000)%>% # 1000 km buffer. 
        st_transform(latlong)%>%
        st_as_text()
      
      bb_arctic <- bioregions%>%
        filter(REGION == "Arctic")%>%
        st_bbox()%>%
        st_as_sfc()%>%
        st_transform(CanProj)%>% #area based projection (units m)
        st_buffer(1000 * 1000)%>% # 1000 km buffer
        st_transform(latlong)%>%
        st_as_text()
      
      bb_pacific <- bioregions%>%
        filter(REGION == "Pacific")%>%
        st_bbox()%>%
        st_as_sfc()%>%
        st_transform(CanProj)%>% #area based projection (units m)
        st_buffer(1000 * 1000)%>% # 1000 km buffer
        st_transform(latlong)%>%
        st_as_text()
      
      bb_zones <- list(bb_atlantic,bb_arctic,bb_pacific)
      names(bb_zones) <- ocean_ord
      
      oceans <- bioregions%>%
        filter(REGION %in% ocean_ord)%>%
        group_by(REGION)%>%
        summarise(geometry=st_union(geometry))%>%
        ungroup()%>%
        st_as_sf()%>%
        st_transform(latlong)%>%
        mutate(region=factor(REGION,levels=ocean_ord))%>% #sort
        arrange(region)%>%
        dplyr::select(-region)
      
      #basemapes
        canada_eez <- read_sf("data/shapefiles/Canada_EEZ.shp")%>%
                      st_transform(CanProj)
        
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
                           mutate(country="Greenland"),
                         ne_states(country = "Russia",returnclass = "sf")%>%
                           dplyr::select(name_en,geometry)%>%
                           st_as_sf()%>%
                           st_union()%>%
                           st_transform(latlong)%>%
                           st_as_sf()%>%
                           mutate(country="Russia"))%>%
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

#Check for missing using CARMs
    ## read in Canadian Register of Marine Species (CARMs) checklist download - for all the species. 
    carms_df <- read.csv("data/CaRMS_checklist_20230614_Canada.csv")%>%
                filter(UnacceptReason == "",
                       Class=="Teleostei")%>% #just fish
                dplyr::select(AphiaID,ScientificName,Locality,
                              Kingdom,Phylum,Class,Order,Family,Genus,
                              Subgenus,Species,Subspecies)
    
    missing_carms_id <- setdiff(unique(carms_df$AphiaID),unique(worms_df$AphiaID)) #in CaRMS but not on our list
             
    missing_carms_species <- carms_df%>%
                        filter(AphiaID %in% missing_carms_id,
                               Species !="")%>%
                        distinct(AphiaID,.keep_all = TRUE)%>%
                        pull(ScientificName)
    
    #run classification based on the species names from WORMs
    missing_carms_classification <- list()
    
    for(i in 1:length(missing_carms_species)){
      
      message(paste0("Working on ",missing_carms_species[i]," ",i,"/",length(missing_carms_species)))
      
      missing_carms_classification[[missing_carms_species[i]]] <- worms_classify(missing_carms_species[i])
      
    }
    
    carms_missing_df <- do.call("rbind",missing_carms_classification)
    
    #save interim output
    save(carms_missing_df,file="output/carms_missing_df.RData")
    
#Identify missing using OBIS records

    #Identify taxa in Canada using the area IDs 
      can_taxa_arctic <- checklist(areaid = 33)%>%as_tibble()%>%mutate(ocean="Arctic") #Arctic
      can_taxa_atlantic <- checklist(areaid = 34)%>%as_tibble()%>%mutate(ocean="Atlantic") #Atlantic
      can_taxa_pacific <- checklist(areaid = 35)%>%as_tibble()%>%mutate(ocean="Pacific") #Pacific
      
      can_taxa <- rbind(can_taxa_arctic,can_taxa_atlantic,can_taxa_pacific)

    #save outputs that can be reloaded instead of re-extracting. 
      save(can_taxa,file=paste0("output/obis_checklist_all.RData"))

      #assemble the data.frame
      can_df <- can_taxa%>%
                filter(taxonomicStatus == 'accepted',#keep if it has a valid iD
                       as.logical(is_marine)|as.logical(is_brackish), #keep if it is marine or brackish
                       superclass=="Actinopteri" | subclass=="Teleostei" | class=="Teleostei",
                       !is.na(species))%>% #keep the bony fishes
                dplyr::select(scientificName,taxonID,bold_id,ncbi_id,kingdom,phylum,class,order,family,genus,species)%>%
                rename(AphiaID = taxonID)
      
      #key out the missing species    
  missing_obis_id <- setdiff(unique(can_df $AphiaID),unique(worms_df$AphiaID))%>%
                     setdiff(.,missing_carms_id)
 
  missing_obis_species <- can_taxa%>%filter(taxonID %in% missing_obis_id)%>%pull(scientificName)

  #report back 
    message(paste0(length(missing_obis_id))," species identifed by OBIS were not listed in the reference database or on CaRMS.")
    message(setdiff(unique(c(worms_df$AphiaID,missing_carms_id)),unique(can_df $AphiaID))%>%length()," species in the reference database and CaRMs, but not listed in OBIS within the Canadian EEZ.")

#run identification script
missing_obis_classification <- list()

for(i in 1:length(missing_obis_species)){
  
  message(paste0("Working on ", missing_obis_species[i]," ",i,"/",length(missing_obis_species)))
  
  missing_obis_classification[[missing_obis_species[i]]] <- worms_classify(missing_obis_species[i])
  
}

obis_missing_df <- do.call("rbind",missing_obis_classification)

#save interim output
save(obis_missing_df,file="output/obis_missing_df.RData")

#Geographic search --------------------------------

    ##Key out the oceans for the CaRMs species not on obis 
    carms_obis_missing <- setdiff(missing_carms_id,missing_obis_id)
    
    sp_names <- carms_df%>%
                filter(AphiaID %in% carms_obis_missing)%>%
                distinct(AphiaID,.keep_all = TRUE)%>%
                filter(!Species == "")%>%
                pull(ScientificName)
    
    eez_thresh <- 250 #distance in km for an obis observation to be associated 

    ocean_intersection <- list()

        for(i in sp_names){
          
          message(paste0("working on "),i," ",which(i == sp_names),"/",length(sp_names)) #progress message 
          
          sp_temp <- NULL
          
          for(j in 1:3){
            
            message(paste0(names(bb_zones)[j])," record search ..")  
            
            temp <- robis::occurrence(taxonid=carms_df%>%filter(ScientificName==i)%>%distinct(AphiaID)%>%pull(AphiaID),geometry = bb_zones[j])
            
            if(length(temp) !=0){sp_temp <- rbind(sp_temp,temp%>%
                                                    mutate(ocean = (names(bb_zones)[j]))%>%
                                                    dplyr::select(scientificName,species,
                                                                  decimalLatitude,decimalLongitude,ocean))}
            
            #save interim outputs
            save(sp_temp,file=paste0("output/robis_records/carms_extras/robis_",gsub(" ","_",i),".RData"))
            
            rm(temp)
            
          } #end 'j' loop
          
          if(!is.null(sp_temp)){
            
            #format data as a spatial data.frame
            sf_format <- sp_temp%>%
                        dplyr::select(ocean,decimalLatitude,decimalLongitude)%>%
                        st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
                        mutate(ocean = ocean_ord[as.numeric(st_intersects(.,oceans))])
                      
            #if no points fall within a 'Canadian' ocean check which is the nearest (as the crow would fly)
            if(sum(is.na(sf_format$ocean)) == nrow(sf_format)){
              
              ocean_dist <- sf_format%>%
                mutate(ocean = ocean_ord[st_nearest_feature(.,oceans)])
             
                  #now get the distance from each record to it's (Canadian) nearest ocean polygon
                  for(j in unique(ocean_dist$ocean)){
                    
                    ocean_dist[j,"dist"] <- st_distance(ocean_dist[j,],oceans%>%
                                                          filter(REGION == ocean_dist[j,]%>%pull(ocean)))%>%as.numeric()/1000
                    
                  } # end j loop
              
              #get rid of any duplicate distances. 
              ocean_dist <- ocean_dist%>%distinct(dist,.keep_all=TRUE) 
              
              #if all distances are further than the limit only include the closest one
              if(min(ocean_dist$dist)>eez_thresh){ocean_dist <- ocean_dist%>%
                filter(dist==min(dist,na.rm=T))}else{ocean_dist <- ocean_dist%>%filter(dist<=eez_thresh)}
              
              out <- data.frame(species=i,
                                ocean = ocean_dist%>%pull(ocean)%>%unique()%>%paste(.,collapse="-"),
                                dist = min(ocean_dist$dist)) #closest observation
              
            }else{
              
              out <- data.frame(species=i,
                                ocean=sf_format%>%filter(!is.na(ocean))%>%pull(ocean)%>%unique()%>%paste(.,collapse="-"),
                                dist=0)
              
            }
            
          }#end if(!is.null(sp_temp))
          
          #no observations within 1000 km of Canadian oceans 
          if(is.null(sp_temp)){
            
            out=data.frame(species=i,
                           ocean=NA,
                           dist=NA)
          }
          
          ocean_intersection[[i]] <- out
          
          rm(out)
          
        } #end 'i' loop

    save(ocean_intersection,file="missing_carms_ocean_intersection.RData")
    
    carms_ocean_df <- do.call("rbind",ocean_intersection)%>%
      arrange(species,ocean)%>%
      mutate(input=species)%>% #match for the species whereby the name from CaRMS doesn't match worms
      left_join(.,carms_missing_df,by="input")%>% 
      dplyr::select(-species.x)%>%
      rename(species=species.y)%>%
      distinct(AphiaID,.keep_all=TRUE)
    filter(!is.na(ocean))
    
    carms_null_obs <- carms_ocean_df%>%filter(is.na(ocean)) #these species have no records within the distance threshold of Canada
    
    carms_ocean_df <- carms_ocean_df%>%filter(!is.na(ocean))
    
    #save interim outputs
    save(carms_ocean_df,file="output/carms_ocean_df.Rdata")
    save(carms_null_obs,file="output/carms_null_obs.RData")
    
###Compile the original worms based database, new IDs from CaRMS and from the OBIS EEZ search. 
    can_spec_list <- rbind(worms_df%>%mutate(source="Org_CanFish"),
                           carms_ocean_df%>%mutate(source="CaRMS")%>%select(c(names(worms_df),"source")),
                           obis_missing_df%>%mutate(source="OBIS")%>%select(c(names(worms_df),"source")))
    
    #save interim output
    save(can_spec_list,file="output/can_spec_list.RData")
                      
## Extended search ------------ 
    
    #adjacent american Atlantic and pacific waters
    us_alaska <- checklist(areaid = 265)%>%as_tibble()%>%mutate(ocean="Alaska",areaid = 265)
    us_north_pacific <- checklist(areaid = 274)%>%as_tibble()%>%mutate(ocean="North Pacific",areaid=274)
    us_north_atlantic <- checklist(areaid = 272)%>%as_tibble()%>%mutate(ocean="North Atlantic",areaid=272)

    us_taxa <- rbind(us_alaska,us_north_pacific,us_north_atlantic)

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
    us_missing <- setdiff(us_df%>%pull(AphiaID),unique(c(worms_df%>%pull(AphiaID),
                                                         missing_carms_id,
                                                         missing_obis_id)))
    
    us_sp <- us_df%>%
             filter(AphiaID %in% us_missing)%>%
             pull(scientificName)
    
    areas <- data.frame(region=c("Alaska","North Pacific","North Atlantic"),
                        areaid=c(265,274,272))
    
    for(i in us_sp){
      
      area_id <- 
      
      temp <- robis::occurrence(taxonid=us_df%>%filter(ScientificName==i)%>%distinct(AphiaID)%>%pull(AphiaID),areaid =)
      
      
    }
    

# #create a simplified buffer for the search within set buffer of Canada. 
# buf <- 250 #kms
# 
# can_buff <- canada_eez%>%
#   st_buffer(buf*1000)%>% #buffer in km
#   # st_bbox()%>%
#   # st_as_sfc()%>%
#   # st_as_sf()%>%
#   st_transform(latlong)%>%
#   st_as_sfc()%>%
#   st_simplify()
# 
# plot_dim <- can_buff%>%
#   st_transform(CanProj)%>%
#   st_buffer((buf+100)*1000)%>%
#   st_bbox()
# 
# #map it out to view the search area
# p1 <- ggplot()+
#   geom_sf(data=can_buff%>%st_transform(CanProj),fill=NA)+
#   geom_sf(data=basemap,fill="grey70")+
#   geom_sf(data=canada_eez,fill=NA)+
#   geom_sf(data=Canada,fill="seagreen4")+
#   theme_bw()+
#   coord_sf(expand=0,xlim=plot_dim[c(1,3)],ylim=plot_dim[c(2,4)]);p1
# 
# ggsave(paste0("output/obis_seach_",buf,"km.png"),p1,width=6,height=6,dpi=300,units="in")
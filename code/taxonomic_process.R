### Vizualizations of coverage

#load libraries -----
    library(dplyr)
    library(ggplot2)
    library(taxize)
    library(tidyr)
    library(sf)
    library(robis)
    library(rfishbase)
    library(rgbif)

#functions -------
id_fun <- function(x,tax_order,return="id"){ #extract the lowest id key 
  
  id <- x%>%
    filter(rank %in% tax_order)%>%
    mutate(rank=factor(rank,levels=tax_order))%>%
    arrange(rank)%>%
    slice(n())%>%
    pull(id)%>%
    as.numeric()
  
  rank <- x%>%
    filter(rank %in% tax_order)%>%
    mutate(rank=factor(rank,levels=tax_order))%>%
    arrange(rank)%>%
    slice(n())%>%
    pull(rank)%>%
    as.character()
  
  if(return=="id"){return(id)}
  if(return=="rank"){return(rank)}
  
}

#projections -----------
    latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
    utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0" #works for the Scotian Shelf
    CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    
#polygons ----------
    bioregions <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/DFO_Marine_Bioregions_Clipped_1M_CAEAC_2012_05_31.shp")
    
#load species list --------
    species_list <- read.csv("data/Coad_OBIS_fish_list_Canada.csv")

## run classification on the species list 
    
    #using worms--
    worms_classification <- list()
    
    for(i in 1:nrow(species_list)){
      
      message(paste0("Working on ",species_list[i,"Species"]," ",i,"/",nrow(species_list)))
      
      worms_classification[[species_list[i,"Species"]]] <- worms_classify(species_list[i,"Species"])
      
      }
    
    worms_df <- do.call("rbind",worms_classification)
    
    ## now re-do the analysis using fuzzy matching for the missing names using the worrms API and wm_records_name  
    null_sp <- data.frame(org=worms_df%>%filter(is.na(AphiaID))%>%pull(input),
                          fixed=NA)
    
    for(i in null_sp$org){
      
      null_sp[null_sp$org==i,"fixed"] <- wm_records_taxamatch(name=i,marine_only=FALSE)%>%
        data.frame()%>%pull(valid_name)
      
    }
    
    worms_classification2 <- list()
    
    for(i in 1:nrow(null_sp)){
      
      message(paste0("Working on ",null_sp[i,"org"]," ",i,"/",nrow(null_sp)))
      
      worms_classification2[[null_sp[i,"org"]]] <- worms_classify(null_sp[i,"fixed"])
      
    }
    
    null_df <- do.call("rbind",worms_classification2)%>%
                mutate(input=null_sp$org) #this will match the original data
    
    #now add those entries to the database. 
    worms_df <- worms_df%>%
                filter(!is.na(AphiaID))%>%
                rbind(.,null_df)
    
    #fix one error for environment
    worms_df%>%filter(env=="") # Synaphobranchus pinnatus -- marine animal
    
    worms_df <- worms_df%>%
                mutate(env=ifelse(env=="","Marine",env))

    save(worms_df,file = "output/taxonomy_wide.RData")

#now identify what ocean they were identified in ------------
        
        ocean_ord <- c("Atlantic","Arctic","Pacific")

        #create bounding boxes that can be used by r-obis with a 1000 km boundary.  
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
        
        #for each species check the occurrence records for each region
        
        ocean_intersection <- list()
        
        sp_names <- worms_df%>%
                    filter(env!="Freshwater")%>% #do not search for the freshwater species. 
                    pull(species)
        
        eez_thresh = 250 #distance in km for an obis observation to be associated 
        
        for(i in sp_names){
          
          message(paste0("working on "),i," ",which(i == sp_names),"/",length(sp_names)) #progress message 
          
          sp_temp <- NULL
          
          for(j in 1:3){
          
            message(paste0(names(bb_zones)[j])," record search ..")  
            
            temp <- robis::occurrence(taxonid=worms_df%>%filter(species==i)%>%pull("AphiaID"),geometry = bb_zones[j])
            
            if(length(temp) !=0){sp_temp <- rbind(sp_temp,temp%>%
                                                    mutate(ocean = (names(bb_zones)[j]))%>%
                                                    dplyr::select(scientificName,species,
                                                                  decimalLatitude,decimalLongitude,ocean))}
            
            #save interim outputs
            save(sp_temp,file=paste0("output/robis_records/robis_",gsub(" ","_",i),".RData"))
            
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
            for(j in 1:nrow(ocean_dist)){
              
              ocean_dist[j,"dist"] <- st_distance(ocean_dist[j,],oceans%>%
                                                    filter(REGION == ocean_dist[j,]%>%pull(ocean)))%>%as.numeric()/1000
              
            }
            
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
      
#save the interium outputs
save(ocean_intersection,file="output/ocean_intersection.RData")   

#bind the output together
ocean_df <- do.call("rbind",ocean_intersection)%>%
            arrange(species,ocean)%>%
            left_join(.,worms_df)


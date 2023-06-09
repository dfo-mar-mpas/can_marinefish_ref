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

        #create bounding boxes that can be used by r-obis 
        bb_atlantic <- bioregions%>%
              filter(REGION == "Atlantic")%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_transform(latlong)%>%
              st_as_text()
        
        bb_arctic <- bioregions%>%
              filter(REGION == "Arctic")%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_transform(latlong)%>%
              st_as_text()
        
        bb_pacific <- bioregions%>%
              filter(REGION == "Pacific")%>%
              st_bbox()%>%
              st_as_sfc()%>%
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
        
        sp_names <- worms_df$species
        
        for(i in sp_names){
          
          message(paste0("working on "),i) #progress message 
          
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
        } #end 'i' loop
        
        #narrow down to within the EEZ and within a threshold of the EEZ
        
        eez_thresh = 250 #distance in km for an obis observation to be associated 
        
        robis_names <- dir("output/robis_records/",full.names = TRUE) #these are the species within the bounding boxes 
        
      ocean_intersection <- list()
      null_extracts <- NULL
        
        for (i in robis_names){
          
          message(paste0("Working on ",gsub("output/robis_records/robis_","",i)," ",which(i == robis_names),"/",length(robis_names)))
          
          load(i) #load the extractions from previous loop
          
          if(!is.null(sp_temp)){
         
          sf_format <- sp_temp%>%
                       dplyr::select(ocean,decimalLatitude,decimalLongitude)%>%
                       st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
                       mutate(ocean = ocean_ord[as.numeric(st_intersects(.,oceans))])
        
            #if no points fall within the ocean find out which ocean is the closest and how close
            if(sum(is.na(sf_format$ocean)) == nrow(sf_format)){
              
              ocean_dist <- sf_format%>%
                mutate(ocean = ocean_ord[st_nearest_feature(.,oceans)])
              
              for(j in 1:nrow(ocean_dist)){
              
              ocean_dist[j,"dist"] <- st_distance(ocean_dist[j,],oceans%>%
                                                                 filter(REGION == ocean_dist[j,]%>%pull(ocean)))%>%as.numeric()/1000

              }
              
              ocean_dist <- ocean_dist%>%distinct(dist,.keep_all=TRUE) #get rid of any duplicate distances. 
              
              #if all distances are further than the limit only include the closest one
              if(min(ocean_dist$dist)>eez_thresh){ocean_dist <- ocean_dist%>%filter(dist==min(dist,na.rm=T))}else{ocean_dist <- ocean_dist%>%filter(dist<=eez_thresh)}
              
              out <- data.frame(species=unique(sp_temp$species),
                                ocean = ocean_dist%>%pull(ocean)%>%unique(),
                                dist = min(ocean_dist$dist)) #closest observation
              
            } #ocean distance check
          
            #if there are observations from within
              if(sum(is.na(sf_format$ocean)) != nrow(sf_format)){
              
                  out <- data.frame(species=unique(sp_temp$species),
                                  ocean = sf_format%>%filter(!is.na(ocean))%>%pull(ocean)%>%unique(),
                                  dist=0)
                
              } #0 demarking 'within' while keeping it numeric 
          
          ocean_intersection[[unique(sp_temp$species)]] <- out
          
          } #null check 
          
          if(is.null(sp_temp)){null_extracts <- c(null_extracts,i)}
          
          rm(sp_temp)
          
    }  #ocean 'i' intersection for each obis extraction 
      
#save the interium outputs
save(ocean_intersection,file="output/ocean_intersection.RData")
save(null_extracts,file="output/null_extracts.RData")      

#bind the output together
ocean_df <- do.call("rbind",ocean_intersection)%>%
  arrange(species,ocean)


table(ocean_df$ocean)



#Do re_extractions at a larger extent for the null_extracts ----------

bb_atlantic_lg <- bioregions%>%
                  filter(REGION == "Atlantic")%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_transform(CanProj)%>% #area based projection (units m)
                  st_buffer(1000 * 1000)%>% # 1000 km buffer. 
                  st_transform(latlong)%>%
                  st_as_text()

bb_arctic_lg <- bioregions%>%
                filter(REGION == "Arctic")%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_transform(CanProj)%>% #area based projection (units m)
                st_buffer(1000 * 1000)%>% # 1000 km buffer
                st_transform(latlong)%>%
                st_as_text()

bb_pacific_lg <- bioregions%>%
                  filter(REGION == "Pacific")%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_transform(CanProj)%>% #area based projection (units m)
                  st_buffer(1000 * 1000)%>% # 1000 km buffer
                  st_transform(latlong)%>%
                  st_as_text()

bb_zones_lg <- list(bb_atlantic_lg,bb_arctic_lg,bb_pacific_lg)
names(bb_zones_lg) <- ocean_ord

null_sp_names <- null_extracts%>%
            gsub("output/robis_records/robis_","",.)%>%
            gsub(".RData","",.)%>%
            gsub("_"," ",.)


for(i in null_sp_names){
 
  message(paste0("working on "),i) #progress message 
  
  sp_temp <- NULL
  
  for(j in 1:3){
    
    message(paste0(names(bb_zones_lg)[j])," record search ..")  
    
    temp <- occurrence(i,geometry = bb_zones_lg[j])
    
    if(length(temp) !=0){sp_temp <- rbind(sp_temp,temp%>%
                                            mutate(ocean = (names(bb_zones_lg)[j]))%>%
                                            dplyr::select(scientificName,species,
                                                          decimalLatitude,decimalLongitude,ocean))}
    
    #save interim outputs
    save(sp_temp,file=paste0("output/robis_records/lgextent/robis_",gsub(" ","_",i),".RData"))
    
    rm(temp)
    
  } #end 'j' loop
} #end 'i' loop


eez_thresh_lg = 1000 #distance in km for an obis observation to be associated 

robis_names_lg <- dir("output/robis_records/lgextent/",full.names = TRUE) #these are the species within the bounding boxes 

ocean_intersection_lg <- list()
null_extracts_lg <- NULL

for (i in robis_names_lg){
  
  message(paste0("Working on ",gsub("output/robis_records/lgextent/robis_","",i)," ",which(i == robis_names_lg),"/",length(robis_names_lg)))
  
  load(i) #load the extractions from previous loop
  
  if(!is.null(sp_temp)){
    
    sf_format <- sp_temp%>%
      dplyr::select(ocean,decimalLatitude,decimalLongitude)%>%
      st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
      mutate(ocean = ocean_ord[as.numeric(st_intersects(.,oceans))])
    
    #if no points fall within the ocean find out which ocean is the closest and how close
    if(sum(is.na(sf_format$ocean)) == nrow(sf_format)){
      
      ocean_dist <- sf_format%>%
        mutate(ocean = ocean_ord[st_nearest_feature(.,oceans)])
      
      for(j in 1:nrow(ocean_dist)){
        
        ocean_dist[j,"dist"] <- st_distance(ocean_dist[j,],oceans%>%
                                              filter(REGION == ocean_dist[j,]%>%pull(ocean)))%>%as.numeric()/1000
        
      }
      
      ocean_dist <- ocean_dist%>%distinct(dist,.keep_all=TRUE) #get rid of any duplicate distances. 
      
      #if all distances are further than the limit only include the closest one
      if(min(ocean_dist$dist)>eez_thresh){ocean_dist <- ocean_dist%>%filter(dist==min(dist,na.rm=T))}else{ocean_dist <- ocean_dist%>%filter(dist<=eez_thresh)}
      
      out <- data.frame(species=unique(sp_temp$species),
                        ocean = ocean_dist%>%pull(ocean)%>%unique(),
                        dist = min(ocean_dist$dist)) #closest observation
      
    } #ocean distance check
    
    #if there are observations from within
    if(sum(is.na(sf_format$ocean)) != nrow(sf_format)){
      
      out <- data.frame(species=unique(sp_temp$species),
                        ocean = sf_format%>%filter(!is.na(ocean))%>%pull(ocean)%>%unique(),
                        dist=0)
      
    } #0 demarking 'within' while keeping it numeric 
    
    ocean_intersection_lg[[unique(sp_temp$species)]] <- out
    
  } #null check 
  
  if(is.null(sp_temp)){null_extracts_lg <- c(null_extracts_lg,i)}
  
  rm(sp_temp)
  
}  #ocean 'i' intersection for each obis extraction 
        
#save the interim output
save(ocean_intersection_lg,file="output/ocean_intersection_lg.RData")
       
#Combine the ocean intersections at the localized and large extents
ocean_df <- do.call("rbind",ocean_intersection_lg)%>%
            rbind(.,do.call("rbind",ocean_intersection))%>%
            arrange(species,ocean)%>%
            left_join(.,sp_fishbase)



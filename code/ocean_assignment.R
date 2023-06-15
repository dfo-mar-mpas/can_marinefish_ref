## Ocean_search function 

ocean_assignment <- function(sp_names,eez_thresh=250,maxsearch=1000){

    #sp_names - a vector of scientific names (genus species)  
    #eez_thresh distance (km) (default: 250) for an obis observation to be associated 
    #max_search - max distance (km) to do an extraction from obis (default: 1000)
  
  #load polygons -------------
  bioregions <- read_sf("data/shapefiles/DFO_Marine_Bioregions_Clipped_1M_CAEAC_2012_05_31.shp")
  
  ocean_ord <- c("Atlantic","Arctic","Pacific")
  
  #search extent polygons 
  bb_atlantic <- bioregions%>%
    filter(REGION == "Atlantic")%>%
    st_bbox()%>%
    st_as_sfc()%>%
    st_transform(CanProj)%>% #area based projection (units m)
    st_buffer(max_search * 1000)%>% # 1000 km buffer. 
    st_transform(latlong)%>%
    st_as_text()
  
  bb_arctic <- bioregions%>%
    filter(REGION == "Arctic")%>%
    st_bbox()%>%
    st_as_sfc()%>%
    st_transform(CanProj)%>% #area based projection (units m)
    st_buffer(max_search * 1000)%>% # 1000 km buffer
    st_transform(latlong)%>%
    st_as_text()
  
  bb_pacific <- bioregions%>%
    filter(REGION == "Pacific")%>%
    st_bbox()%>%
    st_as_sfc()%>%
    st_transform(CanProj)%>% #area based projection (units m)
    st_buffer(max_search * 1000)%>% # 1000 km buffer
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

} #end function 
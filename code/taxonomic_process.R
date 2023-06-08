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

    #tax order
    tax_order <- c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","superclass","class",      
                   "order","suborder","family","subfamily","genus","species")   

    # #itis classificaiton 
    # out <-  classification(species_list$Species,db="itis")
    load("output/robis_classification_itis.RData")
    
    #rbind it together
    out_df1 <- do.call(rbind,out)%>%
      mutate(sp = row.names(.),
             sp = gsub('[[:digit:]]+', '', sp),
             sp = gsub("[[:punct:]]", "", sp))

    #key out species not identified in itis by taxize::classificaiton
    null_sp <- out_df1%>%
      filter(is.na(rank))%>%
      pull(sp)

    #gbif on the unclassified------------
    #out2 <- classification(null_sp,db="gbif")
    load("output/robis_classification_gbif_missing.RData")
    
    out2_df1 <- do.call(rbind,out2)%>%
      mutate(sp = row.names(.),
             sp = gsub('[[:digit:]]+', '', sp),
             sp = gsub("[[:punct:]]", "", sp))
    
    #gbif doesn't return all the variables that itis does. So we have to pad the dataframe. 
    
    gbif_missing <- setdiff(tax_order,unique(out2_df1$rank))
    itis_missing <- setdiff(unique(out2_df1$rank),tax_order)
    
    out2_df2 <- out2_df1%>%
                rbind(.,
                      data.frame(sp=rep(unique(out2_df1$sp),each=length(gbif_missing)),
                                 rank=rep(gbif_missing,length(unique(out2_df1))),
                                 name=NA,
                                 id=NA)%>%dplyr::select(names(out2_df1)))%>%
                filter(rank != itis_missing)%>%
                mutate(rank=factor(rank,levels=tax_order))%>%
                arrange(sp,rank)
    
    #clean rownames
    rownames(out_df1) <- NULL
    rownames(out2_df2) <- NULL

    # ### save interim outputs
    #   save(out,file="output/robis_classification_itis.RData") #save the interim outputs so that you don't have to re-run. 
    #   save(out2,file="output/robis_classification_gbif_missing.RData")

  #identify the keys for each species ---- 
  
  names <- out_df1%>%
           filter(!sp %in% null_sp)%>%
           pull(sp)%>%
           unique()
  
    ##ITIS
      #key out the minimum taxonomic level and the associated id
        out_id <- data.frame(sp=names,id=NA,rank=NA)
      
      #extract IDs ** couldn't get this to work quickly wiht group_by for some reason so this is slow but it works. 
        for(i in 1:nrow(out_id)){out_id[i,c("id")] <- id_fun(out_df1%>%filter(sp==out_id[i,"sp"]),tax_order = tax_order,return="id")
                                 out_id[i,c("rank")] <- id_fun(out_df1%>%filter(sp==out_id[i,"sp"]),tax_order = tax_order,return="rank")}
      
      ##GBIF
        #key out the minimum taxonomic level and the associated id
        out_id_gbif <- data.frame(sp=null_sp,id=NA,rank=NA)
        
        #extract IDs ** couldn't get this to work quickly wiht group_by for some reason so this is slow but it works. 
        for(i in 1:nrow(out_id_gbif)){out_id_gbif[i,c("id")] <- id_fun(out2_df1%>%filter(sp==out_id_gbif[i,"sp"]),tax_order = tax_order,return="id")
                                      out_id_gbif[i,c("rank")] <- id_fun(out2_df1%>%filter(sp==out_id_gbif[i,"sp"]),tax_order = tax_order,return="rank")}
        
        
#assemble dataframe in wide format ----
out_df <- out_df1%>%
          filter(rank %in% tax_order)%>%
          mutate(ord=factor(rank,levels=tax_order))%>%
          arrange(sp,ord)%>%
          dplyr::select(sp,name,rank)%>%
          spread(.,key=rank,value=name)%>%
          dplyr::select(sp,all_of(tax_order))%>%
          left_join(.,out_id%>%dplyr::select(sp,id))%>%
          mutate(db="itis")%>%
          rbind(.,
                out2_df2%>%
                  filter(rank %in% tax_order)%>%
                  mutate(ord=factor(rank,levels=tax_order))%>%
                  arrange(sp,ord)%>%
                  dplyr::select(sp,name,rank)%>%
                  spread(.,key=rank,value=name)%>%
                  dplyr::select(sp,all_of(tax_order))%>%
                  left_join(.,out_id_gbif%>%dplyr::select(sp,id))%>%
                  mutate(db="gbif"))
        
## there are three sub-species of two species groups Oncorhynchus clarkii and Esox americanus so we need the sub-species tag
        out_df$subspecies <- NA
        out_df[out_df$species=="Esox americanus","subspecies"] <- c("Esox americanus americanus","Esox americanus vermiculatus")
        out_df[out_df$species=="Oncorhynchus clarkii","subspecies"] <- c("Oncorhynchus clarkii bouvieri","Oncorhynchus clarkii clarkii","Oncorhynchus clarkii lewisi")

## View the mismatches
        out_df%>%
          mutate(spp = case_when(is.na(subspecies)~species,
                                 !is.na(subspecies)~subspecies),
                 match=sp==spp)%>%
          filter(!match)%>%
          dplyr::select(sp,spp) #note these are spelling errors
        
    #save the interim output
        save(out_df,file = "output/taxonomy_wide.RData")
        
# key out whether the fish is freshwater or marine using rfishbase ------------
        sp_names <- out_df%>%distinct(species,.keep_all=TRUE)%>%pull(species)

        sp_fishbase <- data.frame(species=sp_names,domain=NA,accepted=NA)
        
        for(i in 1:length(sp_names)){
          
          #progress message. 
          message(paste0("Working on ",sp_names[i]," ",i," of ",length(sp_names)))
          
          temp <- rfishbase::species(sp_names[i])%>%
                  suppressMessages() # rfishbase returns - Joining with `by = join_by(SpecCode)`
          
            if(nrow(temp) == 0){sp_fishbase[i,"domain"] = "unspecified"}else{
              sp_fishbase[i,"domain"]  <-  temp%>%
                                           data.frame()%>%
                                           mutate(domain = ifelse(Saltwater==1 | Brack==1,"Marine","Freshwater"))%>% #is the fish within the marine environment. 
                                           pull(domain)
            }
          
        }
        
        #save interim output
        save(sp_fishbase,file="output/sp_fishbase.RData")
        
        #get the ones not in fishbase with exact matches. 
        null_fishbase <- sp_fishbase%>%
                         filter(domain=="unspecified")
        
        
        null_fishbase_habitat <- data.frame(species=null_fishbase$species,domain=NA)
        
        for(i in 1:nrow(null_fishbase)){
          
          message(paste0("working on ",null_fishbase[i,"species"]," ",i," of ",nrow(null_fishbase)))
          
          null_fishbase_habitat[i,c("domain","accepted")] <- gbif_alt_habitat(null_fishbase[i,"species"])[,c("habitat","accepted")]
          
        }
        
        save(null_fishbase_habitat,file="output/null_fishbase_habitat.RData") #you can see which do not have accepted names according to GBIF && fishbase. 
        
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
                  filter(REGION %in% ocean_ord)%
                  group_by(REGION)%>%
                  summarise(geometry=st_union(geometry))%>%
                  ungroup()%>%
                  st_as_sf()%>%
                  st_transform(latlong)%>%
                  mutate(region=factor(REGION,levels=ocean_ord))%>% #sort
                  arrange(region)%>%
                  dplyr::select(-region)
        
        #for each species check the occurrence records for each region
        
        sp_names <- out_df$species
        
        for(i in sp_names){
          message(paste0("working on "),i) #progress message 
          
          sp_temp <- NULL
          
          for(j in 1:3){
          
            message(paste0(names(bb_zones)[j])," record search ..")  
            
            temp <- occurrence(i,geometry = bb_zones[j])
            
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
        
        eez_thresh = 200 #distance in km for an obis observation to be associated 
        
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



gbif_alt_habitat <- function(x){
  
  require(rgbif)
  require(dplyr)
  
  accepted <- TRUE
  
  temp <- name_lookup(x,status="ACCEPTED")[["data"]]%>%
          filter(!is.na(habitats))%>%
          pull(habitats)%>%
          unique()
  
  if(is.logical(temp)){
    
    message("name not accepted")
    
    accepted <- FALSE
    
    temp <- name_lookup(x)[["data"]]%>%
      filter(!is.na(habitats))%>%
      pull(habitats)%>%
      unique() 
    
  }
  
  if(length(temp)>1){
    
    hab_list=paste0(temp,collapse="-")
    
    if(!grepl("marine",tolower(hab_list))){temp <- "Freshwater"} else {temp <- "Marine"}
    
  } else {if(!grepl("marine",tolower(temp))){temp <- "Freshwater"} else {temp <- "Marine"}}
  
  return(data.frame(name=x,habitat=temp,accepted=accepted))
  
  
}
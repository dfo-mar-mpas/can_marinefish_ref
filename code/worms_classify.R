worms_classify <- function(x){
  
  require(worrms)
  require(dplyr)
  require(tidyr)
  
  PhyloNames <- c("kingdom","phylum","class","order","family","genus","species")    
  
  #Worms or worrms seems to be full of errors in situations where the species is there but it doesn't come up
  # so we will first start a test function

  test <- try(wm_records_name(x,marine_only = FALSE),silent=TRUE)%>%class() == "try-error"
  
  #Situations where there is an error and nothing returned from worms. 
  if(sum(test)>0){
    
    message("Error with ",x," no matches found.")
    temp2 <- data.frame(matrix(NA,nrow=1,ncol=length(PhyloNames)+3))
    names(temp2) <- c(PhyloNames,"AphiaID","env","input")
    temp2$input <- x
  }
  
  #Test to for situations where there isn't an error but there is no 'accepted' taxonomy
  if(sum(test)==0){
    
    temp <- worrms::wm_records_name(x,marine_only = FALSE)
    
    test2 <- "accepted" %in% temp$status
    
    if(!test2){
      
      temp <- temp%>%slice(1)
      
      environment_meta <- temp%>%
        dplyr::select(isMarine,isFreshwater,isBrackish)%>%
        gather()%>%
        mutate(value=as.logical(value),
               key=substring(key,3,nchar(key)))%>%
        filter(value)%>%
        pull(key)%>%
        paste(.,collapse="-")
      
      temp2 <- wm_classification(temp$valid_AphiaID)%>%
        data.frame()%>%
        mutate(rank=tolower(rank),
               AphiaID=temp$valid_AphiaID,
               env=environment_meta,
               input=x)%>%
        filter(rank %in% PhyloNames)%>%
        spread(key=rank,value=scientificname)%>%
        dplyr::select(all_of(PhyloNames),AphiaID,env,input)
      
      }
  }
  
  #situations where the name works in worms and there is an 'accepted' nomenclature. 
  if(sum(test)==0 & test2){
  
  temp <- worrms::wm_records_name(x,marine_only = FALSE)
  
  if(sum(temp$status == "accepted")>1){
    message(">1 valid names extracting ",temp%>%filter(status=="accepted")%>%slice(1)%>%pull(valid_name))
    temp <- temp%>%filter(status=="accepted")%>%slice(1)
    }
  
  if(sum(temp$status == "accepted")>0){
    
    temp <- temp%>%filter(status=="accepted")
    
    environment_meta <- temp%>%
                        dplyr::select(isMarine,isFreshwater,isBrackish)%>%
                        gather()%>%
                        mutate(value=as.logical(value),
                               key=substring(key,3,nchar(key)))%>%
                        filter(value)%>%
                        pull(key)%>%
                        paste(.,collapse="-")
    
    temp2 <- wm_classification(temp$AphiaID)%>%
             data.frame()%>%
             mutate(rank=tolower(rank),
                    AphiaID=temp$AphiaID,
                    env=environment_meta,
                    input=x)%>%
             filter(rank %in% PhyloNames)%>%
             spread(key=rank,value=scientificname)%>%
             dplyr::select(all_of(PhyloNames),AphiaID,env,input)
  }
  
  if(sum(temp$status == "accepted")==0){
    
    environment_meta <- temp%>%
      dplyr::select(isMarine,isFreshwater,isBrackish)%>%
      gather()%>%
      mutate(value=as.logical(value),
             key=substring(key,3,nchar(key)))%>%
      filter(value)%>%
      pull(key)%>%
      paste(.,collapse="-")
    
    temp2 <- wm_classification(temp$valid_AphiaID)%>%
              data.frame()%>%
              mutate(rank=tolower(rank),
                     AphiaID=temp$valid_AphiaID,
                     env=environment_meta,
                     input=x)%>%
              filter(rank %in% PhyloNames)%>%
              spread(key=rank,value=scientificname)%>%
              dplyr::select(all_of(PhyloNames),AphiaID,env,input)
    
  }
  
  }#end of the working test
  
  return(temp2)

}

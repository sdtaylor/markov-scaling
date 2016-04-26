library(dplyr)
library(tidyr)

#######################################
#Config

#These are the number of points in transects to average together at each scale. Maxing out at 30 over 11 transects
spatial_scales=c(2,4,8,10,15,30)


#######################################
spp_codes=read.csv('./data/SpeciesCodes.txt', sep='\t') %>%
  select(CODE, AGGREGATED.GROUP)

transect_data=read.csv('./data/Tatoosh_Intertidal_Transitions_Transects.txt', sep='\t') %>%
  select(-X) %>%
  gather(date, species, -Transect, -Point.Number) %>%
  left_join(spp_codes, by=c('species'='CODE')) %>% #Convert species to summarized ones defined in spp file
  select(-species) %>%
  rename(species=AGGREGATED.GROUP) %>%
  filter(species!='X', !is.na(species)) %>%
  mutate(month=substring(date, 1,3), year=as.integer(substring(date, 5,6))) %>% #Make dates usuable
  mutate(year=ifelse(year<90, year+2000,year+1900))


#######################################
#Setup the spatial scales to model at

set_id=1
model_sets=data.frame()
for(this_spatial_scale in spatial_scales){
  num_replicates = floor(30/this_spatial_scale)
  
  replicate_identifier=c() #rep is short for replicate
  for(i in 1:this_spatial_scale){
    replicate_identifier=c(replicate_identifier, 1:num_replicates)
  }
  replicate_identifier = sort(replicate_identifier)
  
  #Fill in when the spatial scale isn't a divisor to 30
  while(length(replicate_identifier)<30){
    replicate_identifier = c(replicate_identifier, -1)
  }
  
  for(this_transect in unique(transect_data$Transect)){
    temp_df=data.frame('Transect'=this_transect, Point.Number=1:30, replicate=replicate_identifier,
                       'spatial_scale'=this_spatial_scale, set=set_id)
    model_sets = model_sets %>%
      bind_rows(temp_df)
  }
  
  set_id=set_id+1
}
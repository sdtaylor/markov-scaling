library(dplyr)
library(tidyr)

#######################################
#Config

#These are the number of points in transects to average together at each scale. Maxing out at 30 over 11 transects
spatial_scales=c(2,4,8,10,15,30)


#######################################
#Load data
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


quad_data=read.csv('./data/Tatoosh_Intertidal_Transitions_Quadrats.txt', sep='\t') %>%
  select(-X) %>%
  gather(date, species, -Quadrat, -Point.Number, -Point.Letter) %>%
  left_join(spp_codes, by=c('species'='CODE')) %>% #Convert species to summarized ones defined in spp file
  select(-species) %>%
  rename(species=AGGREGATED.GROUP) %>%
  filter(species!='X', !is.na(species)) %>%
  mutate(month=substring(date, 1,3), year=as.integer(substring(date, 5,6))) %>% #make dates usable
  mutate(year=ifelse(year<90, year+2000,year+1900))  %>%
  select(-date)

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

rm(this_spatial_scale, num_replicates, replicate_identifier, this_transect, temp_df,i)
#####################################
#Build the matrix model out of the quad data.

#Only doing 1 annual transition instead of the summer/winter transition
quad_data = quad_data %>%
  filter(month %in% c('May','Jul','Jun')) %>%
  select(-month)

quad_data_next_year = quad_data %>%
  mutate(year=year-1) %>%
  rename(species_next_year=species)
  
transitions_temp = quad_data %>%
  left_join(quad_data_next_year, by=c('Quadrat','Point.Number','Point.Letter','year')) %>%
  filter(!is.na(species_next_year)) #Some years are missing. Plus the final year can't be included. 

#Get probabilites for each species to each other species
transitions = data.frame()
for(this_species in unique(transitions_temp$species)){
  x=transitions_temp %>%
    filter(species==this_species) %>%
    count(species_next_year) %>%
    mutate(n=n/sum(n), species=this_species) 
  
  transitions = transitions %>%
    bind_rows(x)
}

#Make it an actual matrix. cols are speces at t, rows are prob. of becoming that speces at t+1
transitions= transitions %>%
  spread(species, n) 

rownames(transitions) = transitions$species_next_year
transitions = transitions %>%
  select(-species_next_year) %>%
  as.matrix()

rm(quad_data_next_year, transitions_temp, this_species, x)
########################################################
#Run the markov model given a set of inital conditions and timesteps.
#Return a timeseries of community composition

run_model=function(model, initial_conditions, timesteps){
  results=matrix(nrow=14, ncol=timesteps)
  community_state=initial_conditions
  for(i in seq(timesteps)){
    community_state = model %*% community_state
    results[,i]=community_state
  }
  colnames(results) = seq(timesteps)
  return(results)
}

######################################################
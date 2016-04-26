library(dplyr)
library(tidyr)
library(magrittr)

#######################################
#Config

#These are the number of points in transects to average together at each scale. Maxing out at 30 over 11 transects
spatial_scales=c(2,4,8,10,15,30)

results_file='./results/results.csv'
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
  mutate(year=ifelse(year<90, year+2000,year+1900)) %>%
  filter(month %in% c('May','Jul','Jun'), year<=2012) %>%
  select(-month, -date) %>%
  mutate(point_id=paste(Transect, Point.Number, sep='-')) %>% #Interpret multiple transects as one giant transect
  select(-Transect, -Point.Number)


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

#list of species needed to fill in various things
all_species=quad_data %>% 
  select(species) %>%
  distinct() %>%
  arrange() %>%
  extract2('species')
#######################################
#Setup the spatial scales to model at

set_id=1
model_sets=data.frame()
all_point_ids=sort(unique(transect_data$point_id))
num_points=length(all_point_ids)

for(this_spatial_scale in spatial_scales){
  num_replicates = floor(num_points/this_spatial_scale)
  
  replicate_identifier=c() #rep is short for replicate
  for(i in 1:this_spatial_scale){
    replicate_identifier=c(replicate_identifier, 1:num_replicates)
  }
  replicate_identifier = sort(replicate_identifier)
  
  #Fill in when the spatial scale isn't a divisor to 30
  while(length(replicate_identifier)<num_points){
    replicate_identifier = c(replicate_identifier, -1)
  }
  

  temp_df=data.frame(point_id=all_point_ids, replicate=replicate_identifier,
                     'spatial_scale'=this_spatial_scale, set=set_id)
  model_sets = model_sets %>%
      bind_rows(temp_df)
  
  
  set_id=set_id+1
}

model_sets = model_sets %>%
  filter(replicate >0)

rm(this_spatial_scale, num_replicates, replicate_identifier, num_points, all_point_ids, temp_df,i)
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
#Metrics for composition comparison

#mean squared error
mse=function(actual, predicted){
  mean((actual-predicted)^2)
}

#Euclidean distance
eucl_dist=function(actual, predicted){
  sqrt(sum((actual-predicted)^2))
}

######################################################
#Compare precicted and observed timeseries' of composition
#returns a dataframe of different accuracies over the time series
compare_composition=function(observed, predicted){
  if(!all(dim(observed)==dim(predicted))){stop('observed and predicted rows and/or columns do not match')}
  
  results=data.frame()
  for(this_col in seq(ncol(observed))){
    this_mse=mse(observed[,this_col], predicted[,this_col])
    this_eucl=eucl_dist(observed[,this_col], predicted[,this_col])
    
    results = results %>%
      bind_rows(data.frame(timestep=this_col, mse=this_mse, eucl_dist=this_eucl))
  }
  
  return(results)
}

#####################################################
#Create a matrix time series of actual compostion for a community given a dataframe
#subset of transect_data

create_composition_timeseries=function(df){
  #A list of years that includes missing ones
  years = min(df$year):max(df$year)
  
  yearly_composition=matrix(nrow=length(all_species), ncol=length(years))
  colnames(yearly_composition) = years

  for(i in seq_along(years)){
    this_year=years[i]
    x=df %>%
      filter(year==this_year)
    
    #Put in NA's for missing years
    if(nrow(x)==0){
      x=rep(NA, length(all_species))
    } else {
      x = x %>%
        group_by(species) %>%
        summarize(n=n()) %>%
        ungroup() %>%
        mutate(percent_cover=n/sum(n)) %>% #Percent cover of species actually present
        right_join(data.frame(species=all_species), by='species') %>% #Add in the rest of the species and fill in 0's for cover
        mutate(percent_cover = ifelse(is.na(percent_cover), 0, percent_cover)) %>%
        arrange(species) %>%
        extract2('percent_cover')
    }
    
    yearly_composition[,i] = x
  }
  
  return(yearly_composition)
}

###########################################################
#Iterate through all model sets. making predictions from the model thru time, comparing results
#and compiling everything in a df
run_analysis=function(){
  final_results=data.frame()
  for(this_set in sort(unique(model_sets$set))){
    this_set_list = model_sets %>%
      filter(set==this_set) 
    
    for(this_replicate in sort(unique(this_set_list$replicate))){

      this_replicate_point_ids= this_set_list %>%
        filter(replicate==this_replicate) %>%
        select(point_id) %>%
        distinct() %>%
        extract2('point_id')
      
      actual_composition=transect_data %>%
        filter(point_id %in% this_replicate_point_ids) %>%
        create_composition_timeseries()
        
      predicted_composition = run_model(transitions, actual_composition[,1], ncol(actual_composition))
      results_this_replicate=compare_composition(actual_composition, predicted_composition)
        
      results_this_replicate$set=this_set
      results_this_replicate$replicate=this_replicate
        
      final_results = final_results %>%
        bind_rows(results_this_replicate)
      
    }
    
  }

  final_results = final_results %>%
    left_join( select(model_sets, set, spatial_scale) %>% distinct(), by='set') %>%
    filter(!is.na(mse)) #some na values from missing year
  
  return(final_results)
}

if(!file.exists(results_file)){
  results=run_analysis()
  write.csv(results, results_file, row.names = FALSE)
} else {
  print('Results file exists. Delete to rerun analysis')
}




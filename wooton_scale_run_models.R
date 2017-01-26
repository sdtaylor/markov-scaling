library(dplyr)
library(tidyr)
library(magrittr)

#######################################
#Config

#These are the number of points in transects to average together at each scale. Maxing out at 30 over 11 transects
#spatial_scales=c(1,5,10)
#spatial_scales=c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70)
spatial_scales=c(1,3,5,10,20,30,40,50,60,70)

#The temporal scales to average average results over.
temporal_scales=c(1,2,5,10)

training_years=1993:2000
testing_year_initial=2002
testing_years=2003:2012

results_file='./results/results_wooton.csv'
#######################################
#Load data. Quadrat data from training_years will be used to fit the models. 
#Transect data from testing_years will be used to verify it. 
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
  filter(month %in% c('May','Jul','Jun'), year<=2012) %>% #Only use summer census data
  select(-month, -date) %>%
  mutate(point_id=paste(Transect, Point.Number, sep='-')) %>% #Interpret multiple transects as one giant transect
  select(-Transect, -Point.Number) %>%
  filter((year %in% testing_years) | (year == testing_year_initial))

#Drop some poorly sampled transect points
transect_data = transect_data %>%
  filter(!point_id %in% c('SFN-18','ES-12','LF-25','G-10'))

quad_data=read.csv('./data/Tatoosh_Intertidal_Transitions_Quadrats.txt', sep='\t') %>%
  select(-X) %>%
  gather(date, species, -Quadrat, -Point.Number, -Point.Letter) %>%
  left_join(spp_codes, by=c('species'='CODE')) %>% #Convert species to summarized ones defined in spp file
  select(-species) %>%
  rename(species=AGGREGATED.GROUP) %>%
  filter(species!='X', !is.na(species)) %>%
  mutate(month=substring(date, 1,3), year=as.integer(substring(date, 5,6))) %>% #make dates usable
  mutate(year=ifelse(year<90, year+2000,year+1900))  %>%
  filter(month %in% c('May','Jul','Jun')) %>%
  select(-date, -month) %>%
  mutate(point_id=paste(Quadrat, Point.Letter, Point.Number, sep='-')) %>% #interpret all points independently 
  select(-Quadrat, -Point.Letter, -Point.Number) %>%
  filter(year %in% training_years)

#list of species needed to fill in various things
all_species=quad_data %>% 
  select(species) %>%
  distinct() %>%
  arrange() %>%
  extract2('species')
all_species = sort(as.character(all_species))
#######################################
#Setup the spatial scales to model at
#A set is a single combination of spatial and temporal grain
#Within a set are replicates of which each point is assigned to.
#At the smallest grain the number of points = the  number of replicates. 
#At the largest grain there is a single replicate which all points as aggregated into.
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
  
  #Repeat for every temporal scale
  for(this_temporal_scale in temporal_scales){
    
    temp_df=data.frame(point_id=all_point_ids, spatial_replicate=replicate_identifier,
                       'spatial_scale'=this_spatial_scale, 'temporal_scale'=this_temporal_scale,
                       set=set_id)
    model_sets = model_sets %>%
      bind_rows(temp_df)
    set_id=set_id+1
    
  }
  
}

model_sets = model_sets %>%
  filter(spatial_replicate >0)

#Do the same  thing with the testing years, assigning each year to a temporal replicate for each grain size
temporal_model_sets = data.frame()
num_timesteps = length(testing_years)
for(this_temporal_scale in temporal_scales){
  num_replicates = floor(num_timesteps/this_temporal_scale)  

  replicate_identifier=c()
  for(i in 1:this_temporal_scale){
    replicate_identifier = c(replicate_identifier, 1:num_replicates)
  }
  replicate_identifier = sort(replicate_identifier)
  
  while(length(replicate_identifier)<num_timesteps){
    replicate_identifier = c(replicate_identifier, -99)
  }
  
  temp_df = data.frame(year = testing_years, 
                       temporal_replicate = replicate_identifier,
                       temporal_scale = this_temporal_scale)
  
  temporal_model_sets = temporal_model_sets %>%
    bind_rows(temp_df)
}

rm(this_spatial_scale, this_temporal_scale, num_replicates, replicate_identifier, num_points, all_point_ids, temp_df,i)
#####################################
#Build the markov matrix model 
# training_data is a df of c('point_id','species','year')
# temporal_scale is integer
build_transition_matrix=function(training_data){
  data_next_year = training_data %>%
    mutate(year=year-1) %>%
    rename(species_next_timestep = species)
  
  transitions_temp = training_data %>%
    left_join(data_next_year, by=c('point_id','year')) %>%
    filter(!is.na(species_next_timestep)) #Some years are missing. Plus the final year can't be included. 
  
  #Get probabilites for each species to each other species
  transitions = data.frame()
  for(this_species in unique(transitions_temp$species)){
    x=transitions_temp %>%
      filter(species==this_species) %>%
      count(species_next_timestep) %>%
      mutate(n=n/sum(n), species=this_species) 
    
    transitions = transitions %>%
      bind_rows(x)
  }
  
  #Make it an actual matrix. cols are speces at t, rows are prob. of becoming that speces at t+1
  transitions= transitions %>%
    spread(species, n) 
  
  rownames(transitions) = transitions$species_next_year
  transitions = transitions %>%
    select(-species_next_timestep) %>%
    as.matrix()
  
  #0 probability for transitions that were never observed
  transitions[is.na(transitions)]=0
  
  return(transitions)
}


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

#########################################################
#Use a monte-carlo method for estimating the cover of a single
#point across the timeseries of the testing data.
#Returns, for each timestep in the future, the probability that
#the point will be of each class. 
filler = expand.grid(timestep = 1:(length(testing_years)),
                     species = all_species)

run_single_point_model = function(initial_species, model){
  num_monte_carlo_runs=500
  results = data.frame()
  for(i in 1:num_monte_carlo_runs) {
    present_species = initial_species
    for(timestep in 1:(length(testing_years)-1)){
      next_species_probabilites = model[,present_species]
      present_species = sample(all_species, 1, prob = next_species_probabilites)
      results = results %>%
        bind_rows(data.frame(timestep=timestep,
                             species=present_species,
                             markov_run=i))
    }
  }
  
  results = results %>%
    group_by(timestep) %>%
    count(species) %>%
    mutate(species_cover = n/num_monte_carlo_runs) %>%
    select(-n)
  
  #Fill in low probability species.
  results = results %>%
    right_join(filler, by=c('timestep','species')) %>%
    mutate(species_cover = ifelse(is.na(species_cover), 0, species_cover))
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

#R^2 from a 1:1 line
#log transform because Marks & Muller-Landau 2007: 10.1126/science.1140190 
obs_pred_square=function(actual, predicted){
  actual=log1p(actual)
  predicted=log1p(predicted)
  1 - (sum((actual - predicted) ** 2) / sum((actual - mean(actual)) ** 2))
}

##########################################################
#Scoring function

score_observed_vs_predicted = function(observed, predicted){
  predicted = predicted %>%
    rename(species_cover_predicted = species_cover)
  
  both = observed %>%
    left_join(predicted, by=c('spatial_replicate','temporal_replicate','species')) %>%
    group_by(spatial_replicate, temporal_replicate) %>%
    summarize(r2=obs_pred_square(actual = species_cover, predicted = species_cover_predicted)) %>%
    ungroup()
  
  #return(mean(both$r2))
  return(both)
}

###########################################################
#Upscaling functions, used for both predictions and observations

spatial_aggregate = function(df, set_list){
  spatial_set_list = set_list %>%
    select(point_id, spatial_replicate)
  
  df = df %>%
    left_join(spatial_set_list, by='point_id') %>%
    group_by(spatial_replicate, timestep, species) %>%
    summarise(species_cover = mean(species_cover)) %>%
    ungroup()
  
  return(df)
}

temporal_aggregate = function(df, set_list){
  temporal_set_list = set_list %>%
    select(timestep, temporal_replicate)
  
  df = df %>%
    left_join(temporal_set_list, by='timestep') %>%
    group_by(spatial_replicate, temporal_replicate, species) %>%
    summarise(species_cover = mean(species_cover)) %>%
    ungroup()

  return(df)
}


###########################################################
#Iterate through all model sets. making predictions from the model thru time, comparing results
#and compiling everything in a df
run_analysis=function(){
  
  markov_model=build_transition_matrix(quad_data)
  
  timeseries_prediction = data.frame()
  timeseries_actual = data.frame()
  
  timestep_to_year = data.frame(timestep=1:(length(testing_years)),
                                year=testing_years)
  
  #The initial conditions of each point to start the model with. 
  initial_conditions = transect_data %>%
    filter(year==testing_year_initial)
  #The observations to compar model predictions with
  testing_data = transect_data %>%
    filter(year %in% testing_years)
  
  #For each point in the test dataset, run the monte-carlo simulation using the markov
  #model to get estimates of cover at each timestep in the future.
  #TODO: Parallelization should be put in here if I feel the need
  for(this_point_id in unique(transect_data$point_id)){
    point_data = testing_data %>%
      filter(point_id==this_point_id) %>%
      left_join(timestep_to_year, by='year') %>%
      select(-year)
    initial_species = initial_conditions$species[initial_conditions$point_id==this_point_id]
    point_timeseries_prediction = run_single_point_model(initial_species, markov_model)
    point_timeseries_prediction$point_id = this_point_id
    
    timeseries_prediction = timeseries_prediction %>%
      bind_rows(point_timeseries_prediction)
    
    #Add in zeros to observation data
    point_data = point_data %>%
      mutate(species_cover = 1.0) %>%
      right_join(filler, by=c('timestep','species')) %>%
      mutate(species_cover = ifelse(is.na(species_cover), 0, species_cover))
    point_data$point_id = this_point_id
    
    timeseries_actual = timeseries_actual %>%
      bind_rows(point_data)
  }
  
  #For each combination of spatial and temporal grain sizes, aggregate the results and
  #get an accuracy metrics.
  final_results = data.frame()
  for(this_set in sort(unique(model_sets$set))){
    this_set_list = model_sets %>%
      filter(set==this_set) 
    
    this_spatial_scale = unique(this_set_list$spatial_scale)
    this_temporal_scale = unique(this_set_list$temporal_scale)
    
    #Aggregate in space
    this_set_predictions = spatial_aggregate(timeseries_prediction, this_set_list)
    this_set_observations = spatial_aggregate(timeseries_actual, this_set_list)
    
    #Aggregate in time
    this_temporal_set_list = temporal_model_sets %>%
      filter(temporal_scale == this_temporal_scale) %>%
      left_join(timestep_to_year, by='year')
    this_set_predictions = temporal_aggregate(this_set_predictions, this_temporal_set_list)
    this_set_observations = temporal_aggregate(this_set_observations, this_temporal_set_list)
    
    results_this_set = score_observed_vs_predicted(this_set_observations, this_set_predictions)
    results_this_set$temporal_scale = this_temporal_scale
    results_this_set$spatial_scale = this_spatial_scale
    
    final_results = final_results %>%
      bind_rows(results_this_set)
    
    
  }
  
  
  
  return(final_results)
  
}

if(!file.exists(results_file)){
  results=run_analysis()
  write.csv(results, results_file, row.names = FALSE)
} else {
  print('Results file exists. Delete to rerun analysis')
}




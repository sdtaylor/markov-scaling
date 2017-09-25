library(tidyverse)
library(doParallel)
#######################################
#Config

testing_years=1993:2012
training_years=1993:2000
testing_year_initial=1993

results_file='./results/intertidal_model_output.csv'

num_procs=2
#######################################
#Load data. Quadrat data from training_years will be used to fit the models. 
#Transect data from testing_years will be used to verify it. 
spp_codes=read.csv('./data/SpeciesCodes.txt', sep='\t', stringsAsFactors = FALSE) %>%
  select(CODE, AGGREGATED.GROUP)

transect_data=read.csv('./data/Tatoosh_Intertidal_Transitions_Transects.txt', sep='\t', stringsAsFactors = FALSE) %>%
  select(-X) %>%
  gather(date, species, -Transect, -Point.Number) %>%
  left_join(spp_codes, by=c('species'='CODE')) %>% #Convert species to summarized ones defined in spp file
  select(-species) %>%
  rename(species=AGGREGATED.GROUP) %>%
  filter(species!='X', !is.na(species)) %>%
  mutate(month=substring(date, 1,3), year=as.integer(substring(date, 5,6))) %>% #Make dates usuable
  mutate(year=ifelse(year<90, year+2000,year+1900)) %>%
  filter(month %in% c('May','Jul','Jun')) %>% #Only use summer census data
  select(-month, -date) %>%
  mutate(point_id=paste(Transect, Point.Number, sep='-')) %>% #Interpret multiple transects as one giant transect
  select(-Transect, -Point.Number) 

#Drop some poorly sampled transect points
transect_data = transect_data %>%
  filter(!point_id %in% c('SFN-18','ES-12','LF-25','G-10','SP-22'))

quad_data=read.csv('./data/Tatoosh_Intertidal_Transitions_Quadrats.txt', sep='\t', stringsAsFactors = FALSE) %>%
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
  select(-Quadrat, -Point.Letter, -Point.Number) 

#list of species needed to fill in various things
all_species=quad_data %>% 
  select(species) %>%
  distinct() %>%
  arrange() %>%
  pull('species')
all_species = sort(as.character(all_species))

#The unique upoints in the testing dataset
all_point_ids=sort(unique(transect_data$point_id))


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
    for(timestep in 1:(length(testing_years))){
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
    mutate(species_probability = n/num_monte_carlo_runs) %>%
    select(-n)
  
  #Fill in low probability species.
  results = results %>%
    right_join(filler, by=c('timestep','species')) %>%
    mutate(species_probability = ifelse(is.na(species_probability), 0, species_probability))
}


###################################################################
#Setup parallel processing
registerDoParallel(makeCluster(num_procs))

###########################################################
#Iterate through all model sets. making predictions from the model thru time, comparing results
#and compiling everything in a df
run_analysis=function(){
  
  markov_model=build_transition_matrix(quad_data)
  
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
  timeseries_prediction = foreach(this_point_id = all_point_ids, .combine = bind_rows, .packages=c('dplyr'),
                                  .export = c('run_single_point_model','testing_years','all_species','filler')) %dopar% {
    print(this_point_id)
    initial_species = initial_conditions$species[initial_conditions$point_id==this_point_id]
    if(purrr::is_empty(initial_species)){
      return(data.frame())
    }
    
    point_timeseries_prediction = run_single_point_model(initial_species, markov_model)
    point_timeseries_prediction$point_id = this_point_id
    return(point_timeseries_prediction)
  }
  
  #For each point compile the true observations
  timeseries_actual = data.frame()
  for(this_point_id in unique(transect_data$point_id)){
    point_data = testing_data %>%
      filter(point_id==this_point_id) %>%
      left_join(timestep_to_year, by='year') %>%
      select(-year)

    #Add in zeros to observation data
    point_data = point_data %>%
      mutate(species_actual = 1.0) %>%
      right_join(filler, by=c('timestep','species')) %>%
      mutate(species_actual = ifelse(is.na(species_actual), 0, species_actual))
    point_data$point_id = this_point_id
    
    timeseries_actual = timeseries_actual %>%
      bind_rows(point_data)
  }
  
  #For each combination of spatial and temporal grain sizes, aggregate the results and
  #get an accuracy metrics.
  final_results = timeseries_actual %>%
    left_join(timeseries_prediction, by=c('species','point_id','timestep'))
  
  #Put the year back in
  final_results = final_results %>%
    left_join(timestep_to_year, by='timestep') %>%
    select(-timestep)
  
  return(final_results)
  
}

if(!file.exists(results_file)){
  results=run_analysis()
  write.csv(results, results_file, row.names = FALSE)
} else {
  print('Results file exists. Delete to rerun analysis')
}




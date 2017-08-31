library(tidyverse)
#The temporal scales to average average results over.
temporal_scales = c(1,2,5,10)

training_years = 1993:2000
testing_year_initial = 2002
testing_years = 2003:2012

#Used in the cost lost model as the count of B for when
#some action has to be taken
count_threshold = 18

results_file='./results/results_intertidal.csv'

###################################################
#Prepare a lookup scale of how to aggregate the different years
#at each scale
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
rm(temp_df, i, this_temporal_scale, num_replicates, replicate_identifier)
#######################################################
#Clean data

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
  select(-month, -date) 

transect_data$species = as.character(transect_data$species)

prior_year_observation = transect_data %>%
  mutate(year=year+1) %>%
  rename(prior_year_obs = species)
focal_species_data = transect_data %>%
  mutate(B = ifelse(species == 'B',1,0)) %>%
  left_join(prior_year_observation, by=c('Transect','Point.Number','year')) %>%
  filter(!year %in% c(1993,2002)) %>% #Two years without prior years
  mutate(prior_year_obs = ifelse(is.na(prior_year_obs), 'Unknown',prior_year_obs))

########################################################
#Run and predict model

training_data = focal_species_data %>%
  filter(year %in% training_years) 
testing_data = focal_species_data %>%
  filter(year %in% testing_years)

model = glm(B ~ prior_year_obs, data = training_data, family=binomial)

testing_data$predicted = predict(model, newdata = testing_data, type = 'response')

#Convert probability of seeing focul species at each point to transect count data
predictions = testing_data %>%
  group_by(year, Transect) %>%
  summarize(count = sum(B), count_predicted=sum(predicted))

##########################################################
#Apply temporal scales

scaled_predictions = temporal_model_sets %>%
  left_join(predictions, by=c('year')) %>%
  group_by(temporal_scale, temporal_replicate, Transect) %>%
  summarize(count = mean(count), count_predicted = mean(count_predicted)) %>%
  ungroup() 

scaled_predictions$timestep = with(scaled_predictions, (temporal_scale * temporal_replicate) - temporal_scale/2)


##################################################
#Apply decision threshold for cost/loss value model

#forecast outcomes given binary predictions
outcomes = data.frame(threshold_actual=c(0,1,1,0),
                      threshold_predicted=c(1,1,0,0),
                      result_type=c('fp','tp','fn','tn'))

scaled_predictions = scaled_predictions %>%
  mutate(threshold_actual    = count >= count_threshold,
         threshold_predicted = count_predicted >= count_threshold) %>%
  left_join(outcomes, by=c('threshold_actual','threshold_predicted'))

write.csv(scaled_predictions, results_file, row.names = FALSE)

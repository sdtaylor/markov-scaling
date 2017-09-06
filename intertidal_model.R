library(tidyverse)
#The temporal scales to average average results over.
temporal_scales = c(1,2,3,4)

training_years = 1993:2000
testing_year_initial = 2002
testing_years = 2003:2012
all_years = 1993:2012

#Used in the cost lost model as the count of B for when
#some action has to be taken
count_threshold = 8
count_max = 215
count_min = 185

#focal_species = c('B')
focal_species = c('FLC','FILR','FLR')

scaled_results_file='./results/results_intertidal.csv'
cost_loss_value_file = 'results/intertidal_scaled_results.csv'
###################################################
#Prepare a lookup table of how to aggregate the different years
#at each scale
temporal_model_sets = read.csv('temporal_model_sets.csv')


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
  mutate(point_id=paste(Transect, Point.Number, sep='-')) %>% #interpret all points independently 
  select(-month, -date, -Transect, -Point.Number) %>%
  mutate(data_source = 'transect')

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
  mutate(data_source = 'quadrat')

all_data = quad_data %>%
  bind_rows(transect_data)
all_data$species = as.factor(all_data$species)


#all_data$species = as.character(all_data$species)

prior_year_observation = all_data %>%
  mutate(year=year+1) %>%
  select(prior_year_obs = species, year, point_id)
main_data = all_data %>%
  left_join(prior_year_observation, by=c('point_id','year')) %>%
  filter(!year %in% c(1993,2002,2015)) %>% #Two years without prior years
  filter(!is.na(prior_year_obs))

########################################################

make_regressive_predictions = (initial_data, model, num_years){
  
}


#Run and predict model
library(randomForest)

#main_data$species = as.factor(main_data$species)
#main_data$prior_year_obs = as.factor(main_data$prior_year_obs)

training_data = main_data %>%
  filter(data_source == 'quadrat') 
testing_data = main_data %>%
  filter(data_source == 'transect') %>%
  select(year, point_id, species)
initial_year_data = main_data %>%
  filter(data_source == 'transect', year==1994) %>%
  select(year, point_id, prior_year_obs)

model = randomForest(species ~ prior_year_obs, data = training_data)

#Make predictions in auto.regressive fashion
final_predictions = data.frame()
for(this_year in 1995:2012){
  initial_year_data$predicted = predict(model, newdata = initial_year_data, type='prob') 
  initial_year_data = initial_year_data %>%
    select(-prior_year_obs)
  
  final_predictions = final_predictions %>%
    bind_rows(initial_year_data)

  initial_year_data = initial_year_data %>%
    mutate(year = year+1) %>%
    rename(prior_year_obs = predicted)
}


#Convert probability of seeing focal species at each point to transect count data
#
#########################################
#This needs to be auto-reggressive
#########################################

for()


predictions = testing_data %>%
  group_by(year) %>%
  summarize(count = sum(focal_species_present), count_predicted=sum(predicted)) %>%
  ungroup()

x = predictions %>% gather(count_type, value, -year)
ggplot(x, aes(x=year, y=value, color=count_type, group=count_type)) + geom_line() + geom_point() +
  geom_hline(yintercept = 8) 

##########################################################
#Apply temporal scales

scaled_predictions = temporal_model_sets %>%
  left_join(predictions, by=c('year')) %>%
  group_by(temporal_scale, temporal_cell_id) %>%
  summarize(count = sum(count), count_predicted = sum(count_predicted)) %>%
  ungroup() 

scaled_predictions$timestep = with(scaled_predictions, (temporal_scale * temporal_cell_id) - temporal_scale/2)

write.csv(scaled_predictions, scaled_results_file, row.names = FALSE)

#################################################
##################################################
#Apply decision threshold for cost/loss value model
scaled_predictions_binary = temporal_model_sets %>%
  left_join(predictions, by=c('year')) %>%
  #mutate(observed  = 1*(count < count_max & count > count_min),
  #       predicted = 1*(count_predicted < count_max & count_predicted > count_min)) %>%
  mutate(observed  = 1*(count >= count_threshold),
         predicted = 1*(count_predicted >= count_threshold)) %>%
  group_by(temporal_scale, temporal_cell_id) %>%
  summarize(observed = max(observed), predicted = max(predicted)) %>%
  ungroup() 

#The per year per cell costs
treatment_cost = 10
possible_loss_costs = 10 / seq(0.11, 0.99, 0.01)

#forecast outcomes given binary predictions
outcomes = data.frame(observed=c(0,1,1,0),
                      predicted=c(1,1,0,0),
                      type=c('fp','tp','fn','tn'))

#Gives cost per yr in terms of prescribed treatment vs actual losses
calculate_cost = function(df, treatment_cost, loss_cost, expense_type){
  temporal_scale = unique(df$temporal_scale)
  
  if(expense_type=='perfect'){
    df$predicted = df$observed
  } else if(expense_type=='always') {
    df$predicted = 1
  } else if(expense_type=='never') {
    df$predicted = 0
  } else if(expense_type=='forecast'){
    #No corrections done here
  } else {
    stop('No forecast type')
  }

  df = df %>%
    left_join(outcomes, by=c('observed','predicted'))
  
  total_tp_fp_years = sum(df$type %in% c('fp','tp')) * temporal_scale
  total_fn_years     = sum(df$type == 'fn') * temporal_scale
  total_cost = (treatment_cost * total_tp_fp_years) + (loss_cost * total_fn_years)
  
  total_years = length(unique(df$temporal_cell_id)) * temporal_scale
  print(paste(total_years))
  
  average_per_year_cost = total_cost / total_years
  return(average_per_year_cost)  
}

###############################################################################
#Calculate cost/loss model curve
cost_loss_values=data.frame()
#for(this_transect in unique(scaled_predictions_binary$Transect)){
  for(this_loss_cost in possible_loss_costs){
    smallest_grain_data = scaled_predictions_binary %>%
      filter(temporal_scale == min(temporal_scales))
    
    smallest_grain_perfect = calculate_cost(smallest_grain_data, treatment_cost = treatment_cost, 
                                              loss_cost = this_loss_cost, expense_type='perfect')
    smallest_grain_cost_never = calculate_cost(smallest_grain_data, treatment_cost = treatment_cost, 
                                               loss_cost = this_loss_cost, expense_type='never')
    smallest_grain_cost_maximimum = min(treatment_cost, smallest_grain_cost_never)
    print(smallest_grain_perfect)
    print(smallest_grain_cost_maximimum)
    for(this_temporal_scale in temporal_scales){
      this_scale_data = scaled_predictions_binary %>%
        filter(temporal_scale == this_temporal_scale)
      
      this_scale_expense = calculate_cost(this_scale_data, treatment_cost = treatment_cost, 
                                                 loss_cost = this_loss_cost, expense_type='forecast')
  
      cost_loss_values = cost_loss_values %>%
        bind_rows(data.frame('a' = treatment_cost / this_loss_cost,
                             'expense_max' = smallest_grain_cost_maximimum,
                             'expense_perfect' = smallest_grain_perfect, 
                             'expense_forecast' = this_scale_expense,
                             'temporal_scale' = this_temporal_scale))
    }
    
  }
#}

cost_loss_values$value = with(cost_loss_values, (expense_max - expense_forecast)/(expense_max - expense_perfect))

x = cost_loss_values %>%
  filter(value > 0) %>%
  group_by(a, temporal_scale) %>%
  summarize(value = mean(value))

ggplot(x, aes(x=a, y=value, group=as.factor(temporal_scale), color=as.factor(temporal_scale))) +
  geom_line() + 
  ylim(0,1)

ggplot(cost_loss_values, aes(x=a, y=value, group=as.factor(temporal_scale), color=as.factor(temporal_scale))) +
  geom_line() + 
  ylim(0,1)

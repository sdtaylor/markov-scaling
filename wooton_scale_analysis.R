library(dplyr)
library(ggplot2)
library(tidyr)

results=read.csv('./results/results_wooton.csv') %>%
  filter(timestep!=9)


###################
#Grid plot of timelags x R2
results_summarized = results %>%
  mutate(year_lag=(timestep*temporal_scale) - (temporal_scale/2))  %>%
  group_by(temporal_scale, spatial_scale, year_lag, model_type) %>%
  summarize(mse=mean(mse), eucl_dist=mean(eucl_dist),r2_sd=sd(r2), r2=mean(r2) ) %>%
  ungroup() %>%
  filter(spatial_scale %in% c(1,5,10,30,60))


ggplot(results_summarized, aes(x=year_lag, y=r2, colour=model_type, group=model_type)) + geom_point()+ geom_line()+
  geom_hline(yintercept=0.9) +
  xlab('Years into future') + ylab('R^2') +
  theme_bw() +
  facet_grid(temporal_scale~spatial_scale, labeller=label_both)


########################################################
#graph of errror( r-squared) at a specific point into the future.
#error on the y-axis, temporal_scale on the x, spatial_scale as different lines
#2nd graph same as 1st, but temporal and spatial scales switched


##Only for errors for a prediction ~8 years into the future
lags_to_keep=data.frame(year_lag=c(7.5, 9, 7.5),
                        temporal_scale=c(1,2,5),
                        keep_time_lag='yes')

spatial_scales_to_keep=c(1,5,10,20,40,60)

#Grid plot of timelags x R2
results_summarized = results %>%
  mutate(year_lag=(timestep*temporal_scale) - (temporal_scale/2))  %>%
  group_by(temporal_scale, spatial_scale, year_lag, model_type) %>%
  summarize(r2=mean(r2)) %>%
  ungroup() %>%
  left_join(lags_to_keep, by=c('year_lag','temporal_scale')) %>%
  filter(model_type=='markov', keep_time_lag=='yes', spatial_scale %in% spatial_scales_to_keep)

ggplot(results_summarized, aes(x=temporal_scale, y=r2, group=as.factor(spatial_scale), color=as.factor(spatial_scale)))+
  geom_line(size=2)+
  geom_point(size=3) +
  scale_color_brewer(palette = 'Set2') + 
  theme(panel.grid.major = element_line(colour = "gray38"), 
        panel.background = element_rect(fill = "gray95"), 
        legend.position = "right", legend.direction = "vertical") +
  labs(title = "Prediction accuracy \n of intertidal markov model", 
       x = "Temporal Grain Size (years)", y = "R^2", 
       colour = " Spatial\n Grain Size\n (pin points)") 


ggplot(results_summarized, aes(x=spatial_scale, y=r2, group=as.factor(temporal_scale), color=as.factor(temporal_scale)))+
  geom_line(size=2)+
  geom_point(size=3) +
  scale_color_brewer(palette = 'Set2') + 
  theme(panel.grid.major = element_line(colour = "gray38"), 
        panel.background = element_rect(fill = "gray95"), 
        legend.position = "right", legend.direction = "vertical") +
  labs(title = "Prediction accuracy \n of intertidal markov model", 
       x = "Spatial Grain Size (pin points)", y = "R^2", 
       colour = " Temporal\n Grain Size\n (years)") 


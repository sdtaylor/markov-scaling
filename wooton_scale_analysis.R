library(dplyr)
library(ggplot2)
library(tidyr)

results=read.csv('./results/results_wooton.csv') %>%
  filter(timestep!=9)

max_observations=329

###################
#Grid plot of timelags x R2
results = results %>%
  mutate(year_lag=(timestep*temporal_scale) - (temporal_scale/2))  %>%
  group_by(temporal_scale, spatial_scale, year_lag, model_type) %>%
  summarize(mse=mean(mse), eucl_dist=mean(eucl_dist),r2_sd=sd(r2), r2=mean(r2) ) %>%
  ungroup() %>%
  filter(spatial_scale %in% c(1,5,10,30,60))


ggplot(results, aes(x=year_lag, y=r2, colour=model_type, group=model_type)) + geom_point()+ geom_line()+
  geom_hline(yintercept=0.9) +
  xlab('Years into future') + ylab('R^2') +
  theme_bw() +
  facet_grid(temporal_scale~spatial_scale, labeller=label_both)

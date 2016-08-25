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


################################################################33
#Skill score by incorporating long term average. 
results=read.csv('./results/results_wooton.csv') %>%
  filter(timestep!=9) %>%
  select(spatial_scale, temporal_scale, timestep, model_type, r2, replicate, set) %>%
  spread(model_type, r2) %>%
  mutate(year_lag=(timestep*temporal_scale) - (temporal_scale/2), skill_score=1-((1-markov)/(1-naive))) %>%
  mutate(skill_score=ifelse(skill_score<(-1), -1, skill_score)) %>%
  group_by(spatial_scale, temporal_scale, timestep, year_lag) %>%
  summarize(skill_score=mean(skill_score)) %>%
  ungroup() 

ggplot(results, aes(x=year_lag, y=skill_score)) + geom_point()+ geom_line()+
  xlab('Years into future') + ylab('Skill') +
  theme_bw() +
  facet_grid(temporal_scale~spatial_scale, labeller=label_both) +
  geom_hline(yintercept=0) 


#All other graphs are pointless
#################################
#Grid cells of spatial x temporal scales with accuracy as color gradient
x=results %>%
  filter(spatial_scale %in% c(1,5,10,25,50,70)) %>%
  group_by(temporal_scale, spatial_scale) %>%
  summarize(mse=mean(mse), eucl_dist=mean(eucl_dist),r2_sd=sd(r2), r2=mean(r2) ) %>%
  ungroup()

x$num_predictions_spatial=(max_observations/x$spatial_scale)/max_observations
x$num_predictions_temporal=(19/x$temporal_scale)/19


main_plot=ggplot(x, aes(as.factor(spatial_scale), as.factor(temporal_scale), fill=r2, label=round(r2,2))) + 
  geom_raster() +
  scale_fill_gradient(low='grey100', high='grey40', limits=c(0.3, 1.0)) + geom_text() +
  xlab('Spatial scale (# points)') + ylab('Temporal scale (years)') +
  theme_bw() +
  theme(legend.position='none')

top_bar=ggplot(x %>% select(spatial_scale, num_predictions_spatial) %>% distinct(), aes(as.factor(spatial_scale), num_predictions_spatial))+
  geom_bar(stat='identity') +
  theme_bw()

right_bar=ggplot(x %>% select(temporal_scale, num_predictions_temporal) %>% distinct(), aes(as.factor(temporal_scale), num_predictions_temporal))+
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip()

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

grid.arrange(top_bar, blankPlot, main_plot, right_bar,
             ncol=2, nrow=2, heights=c(2.5,5))

################################################
#accuracy over time with many spatial_scales

x=results %>%
  filter(spatial_scale %in% c(1,5,10,25,50)) %>%
  group_by(timestep, spatial_scale) %>%
  summarize(mse=mean(mse), eucl_dist=mean(eucl_dist),r2_sd=sd(r2), r2=mean(r2) ) %>%
  mutate(eucl_dist = 1-eucl_dist)

ggplot(x, aes(x=timestep, y=r2, colour=as.factor(spatial_scale), group=as.factor(spatial_scale)))+
  geom_line(size=2) + 
  geom_errorbar(aes(ymax=r2+r2_sd, ymin=r2-r2_sd))
  theme(axis.title = element_text(size = 20), 
    plot.title = element_text(size = 18), 
    plot.background = element_rect(fill = "white", 
        linetype = "solid")) +labs(title = "Markov model accuracy over time", 
    x = "years into future", y = "R^2 ") + theme(panel.background = element_rect(fill = NA)) + theme(panel.grid.major = element_line(colour = "black", 
    size = 0.2), panel.grid.minor = element_line(colour = "black", 
    size = 0.2))

################################################
# spatial scale x accuracy + spatial scale x % predictions
x=results %>%
  filter(timestep==10, spatial_scale %in% c(1,5,10,20,30,40,50,60,70)) %>%
  group_by(spatial_scale) %>%
  summarize(mse=1-mean(mse), eucl_dist=1-mean(eucl_dist), r2_sd=sd(r2), r2=mean(r2)) %>%
  ungroup() %>%
  mutate(percent_potential_predictions=floor(max_observations/spatial_scale)/max_observations) %>%
  mutate(usefulnes=percent_potential_predictions * eucl_dist)

#x$mse=scale(x$mse)

ggplot(x, aes(x=spatial_scale)) +
  #geom_line(aes(y=percent_potential_predictions), linetype='dashed')+ geom_point(aes(y=percent_potential_predictions)) +
  geom_line(aes(y=r2)) + geom_point(aes(y=r2)) +
  geom_point(data=results, aes(x=spatial_scale, y=r2)) + stat_binhex()

with(x, plot(percent_potential_predictions~spatial_scale, type='b', ylab='', xlab='Spatial Scale (# of transect points)'))
with(x, points(spatial_scale, r2, type='b', pch=15))
legend('right',c('R^2 +- sd','% predictions made'), pch=c(15,1))
with(x, arrows(spatial_scale, r2-r2_sd, spatial_scale, r2+r2_sd, spatial_scale, length=0.05, angle=90, code=3))

with(x, plot(spatial_scale, r2, type='b', pch=0))



y=results %>%
  filter(timestep %in% c(1,5,10,15,20)) %>%
  group_by(spatial_scale, timestep) %>%
  summarize(mse=mean(mse), eucl_dist=1-mean(eucl_dist)) %>%
  ungroup() %>%
  mutate(percent_potential_predictions=floor(max_observations/spatial_scale)/max_observations) %>%
  mutate(usefulnes=percent_potential_predictions * eucl_dist)

ggplot(y, aes(x=spatial_scale, y=mse, colour=as.factor(timestep), group=as.factor(timestep)))+
  geom_line(size=2)

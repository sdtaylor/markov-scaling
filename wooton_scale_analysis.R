library(dplyr)
library(ggplot2)
library(tidyr)

results=read.csv('./results/results_wooton.csv') 

results_summarized = results %>%
  group_by(temporal_scale, spatial_scale) %>%
  summarise(mse = mean(sum_sq, na.rm=T), sd=sd(sum_sq, na.rm=T), n=n()) 

results_summarized_short = filter(results_summarized, spatial_scale >=1)
results_summarized_short$spatial_name = factor(results_summarized_short$spatial_scale, levels=c(1,5,10,20), labels=c('Spatial Grain: 1',' Spatial Grain: 5','Spatial Grain: 10','Spatial Grain: 20'), ordered = T)
ggplot(results_summarized_short, aes(x=temporal_scale, y=1-mse))+
  geom_point(size=3) +
  geom_line(linetype=2)+
  facet_wrap(~spatial_name, nrow = 1)  + 
  theme(panel.background = element_rect(fill = NA)) + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        panel.background = element_rect(fill = "gray93"),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 20),
        strip.text.x=element_text(size=18),
        strip.text.y=element_text(size=18)) + theme(legend.position = "none") +labs(x = "Temporal Grain", y = "1 - MSE")

results_summarized_short$temporal_name = factor(results_summarized_short$temporal_scale, levels=c(1,2,5,10), labels=c('Temporal Grain: 1',' Temporal Grain: 2','Temporal Grain: 5','Temporal Grain: 10'), ordered = T)
ggplot(results_summarized_short, aes(x=spatial_scale, y=1-mse))+
  geom_point(size=3) +
  geom_line(linetype=2)+
  #ylim(0,1)+
  facet_wrap(~temporal_name, nrow = 1)  + 
  theme(panel.background = element_rect(fill = NA)) + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        panel.background = element_rect(fill = "gray93"),
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 20),
        strip.text.x=element_text(size=18),
        strip.text.y=element_text(size=18)) + theme(legend.position = "none") +labs(x = "Spatial Grain", y = "1 - MSE")

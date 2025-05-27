# Prediction of world plant species diversity
# From RF (random forest) model by Pitcher and Hartig reanalysis
# Critique: https://arxiv.org/pdf/2411.15105

# Need native_myc_latitude_data_2023.RDS to run
# Original figures with full DB
# Can be run with Ext_native_myc_latitude_data_2023.RDS

############################
###### LOAD PACKAGES #######
############################

library(mgcv); 
library(gridExtra); 
library(betareg); 
library(MASS); 
library(lme4); 
library(lmerTest); 
library(lsmeans); 
library(ggeffects); 
library(spdep); 
library(ggplot2); 
library(ncf); 
library(ape); 
library(sjPlot); 
library(gridExtra); 
library(MuMIn);
library(maps); 
library(sf); 
library(car);
library(viridis);
library(tidyverse);
library(randomForest)

############################
######## READ DATA #########
############################

dat <-  readRDS("data/Ext_native_myc_latitude_data_2023.RDS") %>%
# extended data option: 
#dat <-  readRDS("data/native_myc_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area",  "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range",
           "AM","EM","ORC","NM")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  select(-geology) %>%          
  filter(area > 6) %>% # based on paper Patrick shared
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>% 
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), 
         temp = as.vector(scale(temp)), 
         prec = as.vector(scale(log10((prec)+.01))),
         elev_range = as.vector(scale(log10((elev_range)+.01)))) %>% 
  #filter(!entity_class2 == "Non-oceanic") %>% 
  drop_na() %>% 
  filter(sprich < 10000)

############################
######## RUN MODELS ########
##### GET PREDICTIONS ######
############################

dat.oi <- dat %>% filter(entity_class2 == "Oceanic")

dat.ml <- dat %>% filter(entity_class2 == "Mainland")
rf.mod_sprich <- randomForest(sprich ~ latitude +  longitude, data = dat.ml, ntree = 1000) 

# generate world coverage latlon
grid <- tidyr::expand_grid(lat = seq(-90, 90, 5), lon = seq(-180, 180, 5)) %>%
  rename(latitude = lat, longitude = lon) %>%
  mutate(abslatitude = abs(latitude))

# sprich predictions
pred.g <- cbind(grid, fit = predict(rf.mod_sprich, grid)) %>% rename(sprich_exp_rf = fit) 
pred.oi <- cbind(dat.oi, fit = predict(rf.mod_sprich, dat.oi)) %>% rename(sprich_exp_rf = fit) 

############################
#### PLOT DATA GRID RF #####
############################

# read world data
world <- map_data("world")

# write out
#png("figures/ext_ml_rf_pred_map_diversity.jpg", width = 8, height = 5, units = 'in', res = 300)
png("figures/ml_rf_pred_map_diversity.jpg", width = 8, height = 5, units = 'in', res = 300)
ggplot()+
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray30", alpha = 0.5) +
  geom_point(data = pred.g, aes(x = longitude, y = latitude, color = sprich_exp_rf, fill = sprich_exp_rf, size = sprich_exp_rf),  
             pch = 21, stroke = 0.1) +
  geom_jitter()+
  scale_color_viridis(begin = 0.3, end = 1, alpha = 0.75) +
  scale_fill_viridis(begin = 0.3, end = 1, alpha = 0.5) +
  scale_size_continuous(range = c(1, 8)) +
  xlab("") + ylab ("") +
  theme_minimal() +
  coord_sf(ylim = c(-65, 85), xlim = c(-200, 200), expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(color = "Species richness") +
  guides(size = "none") 
dev.off()

############################
##### PLOT DATA OI RF ######
############################

# read world data
world <- map_data("world")

# write out
#png("figures/ext_oi_rf_pred_map_diversity.jpg", width = 8, height = 5, units = 'in', res = 300)
png("figures/oi_rf_pred_map_diversity.jpg", width = 8, height = 5, units = 'in', res = 300)
ggplot()+
  geom_polygon(dat = world, aes(x = long, y = lat, group = group), fill = "gray30", alpha = 0.5) +
  geom_point(dat = pred.oi, aes(x = longitude, y = latitude, color = sprich_exp_rf, fill = sprich_exp_rf, size = sprich_exp_rf),  
             pch = 21, stroke = 0.1) +
  geom_jitter()+
  scale_color_viridis(begin = 0.3, end = 1, alpha = 0.75) +
  scale_fill_viridis(begin = 0.3, end = 1, alpha = 0.5) +
  scale_size_continuous(range = c(1, 8)) +
  xlab("") + ylab ("") +
  theme_minimal() +
  coord_sf(ylim = c(-65, 85), xlim = c(-200, 200), expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(color = "Species richness") +
  guides(size = "none") 
dev.off()

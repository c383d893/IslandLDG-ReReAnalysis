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
library(tidyverse)
library(randomForest)

options(na.action = "na.fail")

packageVersion(c("mgcv")) #‘1.8.41’
packageVersion(c("gridExtra")) # ‘2.3’
packageVersion(c("betareg")) # ‘3.1.4’
packageVersion(c("MASS")) # ‘7.3.58.1’
packageVersion(c("lme4")) # ‘1.1.31’
packageVersion(c("lmerTest")) # ‘3.1.3’
packageVersion(c("lsmeans")) # ‘2.30.0’
packageVersion(c("ggeffects")) # ‘1.1.4’
packageVersion(c("spdep")) # ‘1.2.7’
packageVersion(c("ggplot2")) #'3.4.0'
packageVersion(c("ncf")) # ‘1.3.2’
packageVersion(c("ape")) # ‘5.6.2’
packageVersion(c("sjPlot")) # ‘2.8.12’
packageVersion(c("gridExtra")) # '2.3'
packageVersion(c("MuMIn")) # ‘1.47.1’
packageVersion(c("tidyverse")) #  ‘1.3.2’
packageVersion(c("maps")) # ‘3.4.1’
packageVersion(c("sf")) # ‘1.0.9’
packageVersion(c("car")) # ‘3.1.1’
packageVersion(c("viridis")) # ‘0.6.2’
packageVersion(c("tidyverse")) # ‘1.3.2’

############################
###### LOAD FUNCTIONS ######
############################


R2_function = function( Y, Y_hat) {
  rss <- sum((Y_hat - Y) ^ 2)
  tss <- sum((Y - mean(Y)) ^ 2)
  rsq <- 1 - rss/tss
  return(rsq)
}

#overdispersion function
Check.disp <- function(mod,dat) {
  N <- nrow(dat)
  p <- length(coef(mod))
  E1 <- resid(mod, type = "pearson")
  Dispersion <- sum(E1^2)/ (N-p)
  return(Dispersion)
}

#RAC function
Spat.cor <- function(mod,dat, dist) {
  coords <- cbind(dat$longitude, dat$latitude)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

#RAC function when locations repeat (shift latlon)
Spat.cor.rep <- function(mod,dat, dist) {
  coords <- cbind(dat$longitude, dat$latitude) + matrix(runif(2*nrow(dat), 0, 0.00001), nrow = nrow(dat), ncol = 2)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

############################
######## READ DATA #########
############################

#dat <-  readRDS("data/Ext_native_poll_latitude_data_2023.RDS") %>%
  # extended data option: 
dat <-  readRDS("data/native_poll_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area",  "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range",  
           "poll.b", "poll.ab")) %>%
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
  filter(!entity_class2 == "Non-oceanic") %>%
  drop_na()

dat2 <-  readRDS("data/Ext_native_poll_latitude_data_2023.RDS") %>%
  # extended data option: 
  # dat <-  readRDS("data/native_poll_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist",   "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec",    
           "poll.b","poll.ab")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  select(-geology) %>%                                                              
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>%    
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>%
  filter(area > 6) %>% # based on paper Patrick shared
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), dist = as.vector(scale(log10((dist)+.01))), elev_range = as.vector(scale(log10((elev_range)+.01))),
         temp = as.vector(scale(temp)), prec = as.vector(scale(log10((prec)+.01)))) %>%
  #filter(!entity_class2 == "Non-oceanic") %>%
  drop_na() 

############################
######## ML PREDICT ########
############################

############################
######### MAINLAND #########
############################
dat.ml <- dat %>% filter(entity_class2=="Mainland")

gam.mod_sprich <- gam(sprich ~ s(abslatitude), family=nb(link="log"), data = dat.ml)
rf.mod_sprich <- randomForest(sprich ~ latitude+longitude,  data = dat.ml) 

summary(gam.mod_sprich)

mod <- gam.mod_sprich
new.dat.sprich <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred_sprich <- predict.gam(mod,newdata = new.dat.sprich, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.sprich$abslatitude) 

# check assumptions
gam.check(gam.mod_sprich)
Check.disp(gam.mod_sprich,dat)

#####################################
####### MAINLAND BIOTIC POLL ########
#####################################
r2_RF = data.frame(Mutalism_type = rep(NA, 2), 
                   Latitude = rep(NA, 2), 
                   Longitude = rep(NA, 2), 
                   Full = rep(NA, 2), 
                   GAM = rep(NA, 2))
gam.mod.poll.b <- gam(poll.b ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 

rf.mod.poll.b <- randomForest(poll.b ~ latitude+longitude,  data = dat.ml, ntree = 1000L) 
r2_RF[1,1] = "poll.b"
r2_RF[1,2] = R2_function(dat.ml$poll.b, predict(randomForest(poll.b ~ latitude, data = dat.ml, ntree = 1000L))) 
r2_RF[1,3] = R2_function(dat.ml$poll.b, predict(randomForest(poll.b ~ longitude, data = dat.ml, ntree = 1000L))) 
r2_RF[1,4] = R2_function(dat.ml$poll.b, predict(rf.mod.poll.b)) 
r2_RF[1,5] = R2_function(dat.ml$poll.b, predict(gam.mod.poll.b, type = "response")) 


summary(gam.mod.poll.b)

mod <- gam.mod.poll.b
new.dat.poll.b <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.poll.b <- predict.gam(mod,newdata = new.dat.poll.b, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.poll.b$abslatitude, poll = "poll.b") 

# check assumptions
gam.check(gam.mod.poll.b)
Check.disp(gam.mod.poll.b,dat)

#####################################
####### MAINLAND ABIOTIC POLL #######
#####################################

gam.mod.poll.ab <- gam(poll.ab ~ s(abslatitude) , family=nb(link="log"), data = dat.ml) 

rf.mod.poll.ab <- randomForest(poll.ab ~ latitude+longitude,  data = dat.ml, ntree = 1000L) 
r2_RF[2,1] = "poll.ab"
r2_RF[2,2] = R2_function(dat.ml$poll.ab, predict(randomForest(poll.ab ~ latitude, data = dat.ml, ntree = 1000L))) 
r2_RF[2,3] = R2_function(dat.ml$poll.ab, predict(randomForest(poll.ab ~ longitude, data = dat.ml, ntree = 1000L))) 
r2_RF[2,4] = R2_function(dat.ml$poll.ab, predict(rf.mod.poll.ab)) 
r2_RF[2,5] = R2_function(dat.ml$poll.ab, predict(gam.mod.poll.ab, type = "response")) 

summary(gam.mod.poll.ab)

mod <- gam.mod.poll.ab
new.dat.poll.ab <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.poll.ab <- predict.gam(mod, newdata = new.dat.poll.ab, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.poll.ab$abslatitude, poll = "poll.ab") 

# check assumptions
gam.check(gam.mod.poll.ab)
Check.disp(gam.mod.poll.ab,dat)

############################
########## PLOT ############
##### LAT BY POLL TYPE #####
############################

pred.mainland <- rbind(pred.poll.b,pred.poll.ab)

dat.ml.cond <- dat.ml %>%
  select(abslatitude, poll.b, poll.ab) %>%
  gather(key = "poll", value = "sprich", poll.b, poll.ab)

# create a custom color scale
colScale <- scale_colour_manual(values=c("darkgrey", "deeppink3"))
fillScale <- scale_fill_manual(values=c("darkgrey","deeppink3"))

pred.mainland$poll <- ordered(pred.mainland$poll, levels = c("poll.b", "poll.ab"))
dat.ml.cond$poll <- ordered(dat.ml.cond$poll, levels = c("poll.b", "poll.ab"))

# rename poll for legend:
pred.mainland <- pred.mainland %>% mutate(poll = case_when(poll=="poll.b" ~ "Biotic",
                                                           poll=="poll.ab" ~ "Abiotic"))

dat.ml.cond <- dat.ml.cond %>% mutate(poll = case_when(poll=="poll.b" ~ "Biotic",
                                                       poll=="poll.ab" ~ "Abiotic"))
lat.polltype <-
  ggplot(pred.mainland, aes(x = abslatitude, y = fit,color = poll, fill = poll))+
  geom_line(size = 1) +
  geom_point(data = dat.ml.cond, aes(x = abslatitude, y = sprich, color = factor(poll)), alpha = 0.3, size = 4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,2000)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="Pollination \nSyndrome",override.aes=list(fill =NA,size =3, alpha =0.7,linetype = c(0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/Poll_LatbyPollType.jpg", width=10, height= 10, units='in', res=300)
lat.polltype
dev.off()

############################
###### PREDICT EXP IS ######
####### ALL ISLANDS ########
############################

dat.is.min <- dat %>% filter(entity_class2=="Oceanic") %>% select(c('entity_ID',"abslatitude"))
dat.is.min2 <- dat %>% filter(entity_class2=="Oceanic") %>% select(c('entity_ID',"abslatitude", "latitude", "longitude"))


pred_sprich <- predict.gam(gam.mod_sprich,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_sprich_df <- cbind(dat.is.min,pred_sprich) %>% rename(sprich_exp = fit)
pred_sprich_df_rf <- cbind(dat.is.min,fit=predict(rf.mod_sprich, dat.is.min2)) %>% rename(sprich_exp_rf = fit)


pred_poll.b <- predict.gam(gam.mod.poll.b,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_poll.b_df <- cbind(dat.is.min,pred_poll.b) %>% rename(poll.b_exp = fit)
pred_poll.b_df_rf <- cbind(dat.is.min,fit=predict(rf.mod.poll.b, dat.is.min2)) %>% rename(poll.b_exp_rf = fit)


pred_poll.ab <- predict.gam(gam.mod.poll.ab,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_poll.ab_df <- cbind(dat.is.min,pred_poll.ab) %>% rename(poll.ab_exp = fit)
pred_poll.ab_df_rf <- cbind(dat.is.min,fit=predict(rf.mod.poll.ab, dat.is.min2)) %>% rename(poll.ab_exp_rf = fit)


pred_is.dat <- dat %>%
  filter(entity_class2=="Oceanic") %>%
  left_join(pred_sprich_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_poll.b_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_poll.ab_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_sprich_df_rf, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_poll.b_df_rf, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_poll.ab_df_rf, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('poll.b_exp','poll.ab_exp','sprich_exp', 'poll.b_exp_rf','poll.ab_exp_rf','sprich_exp_rf'), as.integer) %>%
  mutate(propb_exp = poll.b_exp/(poll.b_exp + poll.ab_exp), propb_exp_rf = poll.b_exp_rf/(poll.b_exp_rf + poll.ab_exp_rf)) |> 
  mutate(propb_obs = poll.b/(poll.b + poll.ab))

# write out
saveRDS(pred_is.dat,"data/GAMexp_native_poll_latitude_data_2023.RDS")
saveRDS(r2_RF, "data/r2_RF_Poll.RDS")

############################
####### MAIN MODELS ########
############################

############################
######## CREATE DATA #######
############################

# read data and calc debts
full.dat <- dat2 %>%
  select(c("entity_ID", "dist", "area"))

# version 1; debt & C_debt neg to 0, >1 to 1:
pred_is.dat.shrink <- readRDS("data/GAMexp_native_poll_latitude_data_2023.RDS") %>%
  select(-area) %>%
  mutate(sprich = poll.b + poll.ab, sprich_exp = poll.b_exp + poll.ab_exp) %>%
  mutate(poll.b.diff = poll.b_exp - poll.b, poll.ab.diff = poll.ab_exp - poll.ab, sprichdiff = sprich_exp - sprich) %>%
  mutate(poll.b.debt = (poll.b.diff/poll.b_exp), poll.ab.debt = (poll.ab.diff/poll.ab_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(poll.b.debt = ifelse(poll.b.debt < 0, 0, poll.b.debt)) %>% mutate(poll.ab.debt = ifelse(poll.ab.debt <0, 0, poll.ab.debt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(poll.b.debt = ifelse(poll.b.debt > 1, 1, poll.b.debt)) %>% mutate(poll.ab.debt = ifelse(poll.ab.debt >1, 1, poll.ab.debt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_poll.b.debt = (poll.b.diff/sprichdiff), C_poll.ab.debt = (poll.ab.diff/sprichdiff)) %>% 
  mutate(C_poll.b.debt = ifelse(C_poll.b.debt < 0, 0, C_poll.b.debt)) %>% mutate(C_poll.ab.debt = ifelse(C_poll.ab.debt <0, 0, C_poll.ab.debt)) %>% 
  mutate(C_poll.b.debt = ifelse(C_poll.b.debt > 1, 1, C_poll.b.debt)) %>% mutate(C_poll.ab.debt = ifelse(C_poll.ab.debt >1, 1, C_poll.ab.debt)) %>% 
  left_join(full.dat, by = "entity_ID") 

# version 2; remove negatives: 4% data lost
pred_is.dat.drop <- readRDS("data/GAMexp_native_poll_latitude_data_2023.RDS") %>%
  select(-area) %>%
  mutate(sprich = poll.b + poll.ab, sprich_exp = poll.b_exp + poll.ab_exp) %>%
  mutate(poll.b.diff = poll.b_exp - poll.b, poll.ab.diff = poll.ab_exp - poll.ab, sprichdiff = sprich_exp - sprich) %>%
  mutate(poll.b.debt = (poll.b.diff/poll.b_exp), poll.ab.debt = (poll.ab.diff/poll.ab_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(C_poll.b.debt = (poll.b.diff/sprichdiff), C_poll.ab.debt = (poll.ab.diff/sprichdiff)) %>% 
  filter(poll.b.diff > 0 & poll.ab.diff >0 & sprichdiff > 0) %>%
  filter(poll.b.debt > 0 & poll.ab.debt >0 & Tdebt > 0) %>%
  filter(C_poll.b.debt > 0 & C_poll.ab.debt) %>%
  left_join(full.dat, by = "entity_ID") 

############################
######### PREP DATA ########
############################

# ! choose one
pred_is.dat <- pred_is.dat.shrink
pred_is.dat <- pred_is.dat.drop

pred_is.dat.alld <- pred_is.dat %>% 
  select(c('entity_ID','poll.b.debt','poll.ab.debt')) %>%
  gather(key = "poll", value = "debt", poll.b.debt, poll.ab.debt) %>%
  mutate(poll = case_when(poll == "poll.b.debt" ~ "poll.b",                                   
                             poll == "poll.ab.debt" ~ "poll.ab")) 

pred_is.dat.allcd <- pred_is.dat %>% 
  select(c('entity_ID','C_poll.b.debt','C_poll.ab.debt')) %>%
  gather(key = "poll", value = "debt.c", C_poll.b.debt, C_poll.ab.debt)%>%
  mutate(poll = case_when(poll == "C_poll.b.debt" ~ "poll.b",                                   
                          poll == "C_poll.ab.debt" ~ "poll.ab"))

pred_is.dat.all.diff <- pred_is.dat %>% 
  select(c('entity_ID','sprich','latitude','longitude','abslatitude','poll.b.diff','poll.ab.diff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "poll", value = "diff", poll.b.diff, poll.ab.diff) %>%
  mutate(poll = case_when(poll == "poll.b.diff" ~ "poll.b",                                   
                          poll == "poll.ab.diff" ~ "poll.ab")) 

pred_is.dat.all.exp <- pred_is.dat %>% 
  select(c('entity_ID','poll.b_exp','poll.ab_exp')) %>%
  gather(key = "poll", value = "exp", poll.b_exp, poll.ab_exp) %>%
  mutate(poll = case_when(poll == "poll.b_exp" ~ "poll.b",                                   
                          poll == "poll.ab_exp" ~ "poll.ab")) 

pred_is.dat.all.obs <- pred_is.dat %>% 
  select(c('entity_ID', 'poll.b', 'poll.ab')) %>%
  gather(key = "poll", value = "obs", poll.b, poll.ab) 

pred_is.dat.all.T <- pred_is.dat %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("poll", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(poll = as.factor(poll)) 
  
pred_is.dat.all <- within(pred_is.dat.all, poll <- relevel(poll, ref = "poll.ab"))

############################
######### CUT DATA #########
############################

# set vars
x_var = pred_is.dat.all$abslatitude
poll_var= pred_is.dat.all$poll

y_var = pred_is.dat.all$debt

df <- data.frame("xvar" = x_var,"yvar" = y_var,"poll"=poll_var)

# create binned y-values for the x-axis
quantiles_for_cutting <- quantile(df$xvar,seq(0,1,.20))

# cut the data
df$cuts_raw <- cut(df$xvar,breaks = quantiles_for_cutting, include.lowest = T)

# calculate the average value within each bin of the x-axis data
mean_per_cut <- df %>% group_by(cuts_raw) %>%
  summarize(mean_cut_xval = mean(xvar))

#now use the "mean_per_cut" to define our new cuts
df$cuts_labeled <- as.numeric(as.character(cut(df$xvar,breaks = quantiles_for_cutting,
                                               labels = mean_per_cut$mean_cut_xval, include.lowest = T)))

#now calculate the mean response variable values within each bin: 95% CIs assuming a normal distribution here
aggregated_data <- df %>% group_by(cuts_raw,cuts_labeled,poll) %>%
  summarize(mean_y = mean(yvar),
            sd_y = sd(yvar),
            n_y = length(yvar)) %>%
  mutate(se_y = sd_y/sqrt(n_y),
         low95CI_y = mean_y-1.96*se_y,
         high95CI_y = mean_y+1.96*se_y)

debt.aggregated_data <- aggregated_data

############################
######### ONE MODEL ########
############################

############################
######## WITHIN DEBT #######
############################

pred.allcatd <- glm(debt ~ abslatitude*poll + area + dist +elev_range + prec + (1|entity_ID),weights = debt.weights, data = pred_is.dat.all) 
summary(pred.allcatd)
rac <- Spat.cor.rep(pred.allcatd,pred_is.dat.all,2000)
pred.allcatd.rac  <- glm(debt ~ abslatitude*poll + area + dist +elev_range+ prec + rac +(1|entity_ID),weights = debt.weights, data = pred_is.dat.all) 
summary(pred.allcatd.rac)

pred_is.dat.all$rac <- rac
ref<-lsmeans(pred.allcatd.rac, pairwise ~ poll, data = pred_is.dat.all)

#contrasts
means <- emmeans(pred.allcatd.rac, ~ poll, data = pred_is.dat.all)

#look at means order to determine how to write contrasts:
means

#write contrasts:
contrasts <- list("abiotic v biotic" = c(-1,1)
)

#extract results to df:
results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

# glm diagnostics:
par(mfrow=c(3,2))
# homogenetity of variance:
plot(fitted(pred.allcatd.rac) ~ resid(pred.allcatd.rac, type = "pearson"))
# independence:
plot(pred_is.dat.all$area ~ resid(pred.allcatd.rac, type = "pearson"))
plot(pred_is.dat.all$dist ~ resid(pred.allcatd.rac, type = "pearson"))
plot(pred_is.dat.all$elev_range ~ resid(pred.allcatd.rac, type = "pearson"))
boxplot(resid(pred.allcatd.rac, type = "pearson") ~ pred_is.dat.all$poll)
# do not look for normality:
# outliers:
plot(cooks.distance(pred.allcatd.rac), type ="h")
# check dispersion:
Check.disp(pred.allcatd.rac, pred_is.dat.all)

############################
########## DEBT C ##########
############################

pred.allcatcd <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*poll + area + dist +elev_range+ prec+ (1|entity_ID),weights = debt.c.weights, data = pred_is.dat.all) 
rac <- Spat.cor.rep(pred.allcatcd,pred_is.dat.all,2000)
pred.allcatcd.rac  <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*poll + area + dist +elev_range+ prec+ rac +(1|entity_ID),weights = debt.c.weights, data = pred_is.dat.all) 
summary(pred.allcatcd.rac)

# glm diagnostics:
par(mfrow=c(3,2))
# homogenetity of variance:
plot(fitted(pred.allcatcd.rac) ~ resid(pred.allcatcd.rac, type = "pearson"))
# independence:
plot(pred_is.dat.all$area ~ resid(pred.allcatcd.rac, type = "pearson"))
plot(pred_is.dat.all$dist ~ resid(pred.allcatcd.rac, type = "pearson"))
plot(pred_is.dat.all$elev_range ~ resid(pred.allcatcd.rac, type = "pearson"))
boxplot(resid(pred.allcatcd.rac, type = "pearson") ~ pred_is.dat.all$poll)
# do not look for normality:
# outliers:
plot(cooks.distance(pred.allcatcd.rac), type ="h")
# check dispersion:
Check.disp(pred.allcatcd.rac, pred_is.dat.all)

############################
########### PLOT ###########
##### DEBT:FULL MODEL ######
############################

new.dat.poll.b <- pred_is.dat.all %>% filter(poll=="poll.b") %>%
  mutate(area = mean(pred.allcatd.rac$model$area), dist = mean(pred.allcatd.rac$model$dist),
         elev_range = mean(pred.allcatd.rac$model$elev_range), prec = mean(pred.allcatd.rac$model$prec), rac = mean(pred.allcatd.rac$model$rac))
pred.poll.b <- predict.glm(pred.allcatd.rac,newdata = new.dat.poll.b, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.poll.b$abslatitude) 

new.dat.poll.ab <- pred_is.dat.all %>% filter(poll=="poll.ab") %>%
  mutate(area = mean(pred.allcatd.rac$model$area), dist = mean(pred.allcatd.rac$model$dist),
         elev_range = mean(pred.allcatd.rac$model$elev_range), prec = mean(pred.allcatd.rac$model$prec), rac = mean(pred.allcatd.rac$model$rac))
pred.poll.ab <- predict.glm(pred.allcatd.rac,newdata = new.dat.poll.ab, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.poll.ab$abslatitude) 

allcatd.plot <- 
  ggplot() +
  #B section
  geom_line(data = pred.poll.b, mapping = aes(x = abslatitude, y = fit), color ="deeppink3")+
  geom_ribbon(data = pred.poll.b, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "deeppink3", alpha = 0.5) +
  #geom_point(data = pred_is.dat.all %>% filter(poll=="poll.b"), aes(x = abslatitude, y=debt),color ="deeppink3", alpha= 0.3, size=5) +
  geom_point(data = debt.aggregated_data %>% filter(poll=="poll.b"),aes(x = cuts_labeled-2, y = mean_y),color = "deeppink3",size = 4,shape = 16, alpha = 0.8) +
  geom_errorbar(data = debt.aggregated_data %>% filter(poll=="poll.b"),aes(x = cuts_labeled-2, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "deeppink3", width=4,size=2, alpha = 0.8) +
  #AB section
  geom_line(data = pred.poll.ab, mapping = aes(x = abslatitude, y = fit), color ="darkgrey")+
  geom_ribbon(data = pred.poll.ab, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkgrey", alpha = 0.5) +
  #geom_point(data = pred_is.dat.all %>% filter(poll=="poll.ab"), aes(x = abslatitude, y=debt),color ="darkgrey", alpha= 0.3, size=5) +
  geom_point(data = debt.aggregated_data %>% filter(poll=="poll.ab"),aes(x = cuts_labeled+2, y = mean_y),color = "darkgrey",size = 4,shape = 16, alpha = 0.8) +
  geom_errorbar(data = debt.aggregated_data %>% filter(poll=="poll.ab"),aes(x = cuts_labeled+2, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "darkgrey", width=4,size=2, alpha = 0.8) +
  theme_classic(base_size = 40) +
  #geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylab("Proportional species deficit") +
  xlab("Absolute latitude") +
  ylim(0.4,1.05)

# write out
png("figures/Poll_LatBox_withindebt_fullmodelshrink.jpg", width = 10, height = 10, units = 'in', res = 300)
allcatd.plot
dev.off()

allcatd.plot.points <- 
  ggplot() +
  #B section
  geom_line(data = pred.poll.b, mapping = aes(x = abslatitude, y = fit), color ="deeppink3")+
  geom_ribbon(data = pred.poll.b, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "deeppink3", alpha = 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(poll=="poll.b"), aes(x = abslatitude, y=debt),color ="deeppink3", alpha= 0.3, size=5) +
  #AB section
  geom_line(data = pred.poll.ab, mapping = aes(x = abslatitude, y = fit), color ="darkgrey")+
  geom_ribbon(data = pred.poll.ab, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkgrey", alpha = 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(poll=="poll.ab"), aes(x = abslatitude, y=debt),color ="darkgrey", alpha= 0.3, size=5) +
  theme_classic(base_size = 40) +
  #geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylab("Proportional species deficit") +
  xlab("Absolute latitude") +
  ylim(0,1.05)

# write out
png("figures/Poll_LatPoints_withindebt_fullmodelshrink.jpg", width = 10, height = 10, units = 'in', res = 300)
allcatd.plot.points
dev.off()

############################
########### PLOT ###########
#### C DEBT:FULL MODEL #####
############################

new.dat.poll.b <- pred_is.dat.all %>% filter(poll=="poll.b") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.poll.b <- predict.glm(pred.allcatcd.rac,newdata = new.dat.poll.b, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.poll.b$abslatitude) 

new.dat.poll.ab <- pred_is.dat.all %>% filter(poll=="poll.ab") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.poll.ab <- predict.glm(pred.allcatcd.rac,newdata = new.dat.poll.ab, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.poll.ab$abslatitude) 

allcatcd.plot <- 
  ggplot() +
  #AM section
  geom_line(data = pred.poll.b, mapping = aes(x = abslatitude, y = fit), color ="deeppink3")+
  geom_ribbon(data = pred.poll.b, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "deeppink3", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(poll=="poll.b"), aes(x = abslatitude, y=debt.c),color ="deeppink3", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(poll=="poll.b"),aes(x = cuts_labeled, y = mean_y),color = "deeppink3",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(poll=="poll.b"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "deeppink3", width=3,size=1.5) +
  #EM section
  geom_line(data = pred.poll.ab, mapping = aes(x = abslatitude, y = fit), color ="darkgrey")+
  geom_ribbon(data = pred.poll.ab, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkgrey", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(poll=="poll.ab"), aes(x = abslatitude, y=debt.c),color ="darkgrey", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(poll=="poll.ab"),aes(x = cuts_labeled, y = mean_y),color = "darkgrey",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(poll=="poll.ab"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "darkgrey", width=3,size=1.5) +
  theme_classic(base_size = 15) +
  #geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylab("Contribution deficit") +
  xlab("Absolute latitude") +
  ylim(0,1.2)

# write out
png("figures/Poll_LatPoly_contdebt_fullmodel_shrink.jpg", width = 6, height = 6, units ='in', res = 300)
allcatcd.plot
dev.off()

############################
###### MODEL PER POLL ######
############################

pred_is.dat.poll.b <- pred_is.dat.all %>% filter(poll == "poll.b")
pred_is.dat.poll.ab <- pred_is.dat.all %>% filter(poll == "poll.ab")

############################
########## DEBT C ##########
############################

# poll.b Poly
pred.allcatcd.poll.b <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*area + poly(abslatitude,2,raw = TRUE)*dist +poly(abslatitude,2,raw = TRUE)*elev_range + poly(abslatitude,2,raw = TRUE)*prec, weights = debt.c.weights, data = pred_is.dat.poll.b) 
rac <- Spat.cor.rep(pred.allcatcd.poll.b, pred_is.dat.poll.b, 2000)
pred.allcatcd.poll.b.rac  <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*area + poly(abslatitude,2,raw = TRUE)*dist +poly(abslatitude,2,raw = TRUE)*elev_range + poly(abslatitude,2,raw = TRUE)*prec + rac , weights = debt.c.weights, data = pred_is.dat.poll.b) 
summary(pred.allcatcd.poll.b.rac)
pred.allcatcd.poll.b.rac.min  <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE) + rac , weights = debt.c.weights, data = pred_is.dat.poll.b) 
summary(pred.allcatcd.poll.b.rac.min)

# poll.ab Poly
pred.allcatcd.poll.ab <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*area + poly(abslatitude,2,raw = TRUE)*dist +poly(abslatitude,2,raw = TRUE)*elev_range + poly(abslatitude,2,raw = TRUE)*prec, weights = debt.c.weights, data = pred_is.dat.poll.ab) 
rac <- Spat.cor.rep(pred.allcatcd.poll.ab, pred_is.dat.poll.ab, 2000)
pred.allcatcd.poll.ab.rac  <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*area + poly(abslatitude,2,raw = TRUE)*dist +poly(abslatitude,2,raw = TRUE)*elev_range + poly(abslatitude,2,raw = TRUE)*prec +rac , weights = debt.c.weights, data = pred_is.dat.poll.ab) 
summary(pred.allcatcd.poll.ab.rac)
pred.allcatcd.poll.ab.rac.min  <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE) + rac , weights = debt.c.weights, data = pred_is.dat.poll.ab) 
summary(pred.allcatcd.poll.ab.rac.min)
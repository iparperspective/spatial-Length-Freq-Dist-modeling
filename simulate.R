


########### rspde
library(INLA)
library(inlabru)
library(tidyverse)
source("~/spde-book-functions.R")

rho_size <- 0.75 # size cor
range = 50
params <- c(variance = .5, kappa = 1)


dims = 101
coordin = expand.grid(1:dims,1:dims);names(coordin)=c("x","y");coordin=as.matrix(coordin)
mat = matrix(1,ncol=dims,nrow=dims)
n_size = 16
n = dims^2

max.edge = dims*.05
bound.outer = .3
mesh = inla.mesh.2d(loc.domain = cbind(c(0,dims,dims,0,0),c(0,0,dims,dims,0)),
                    #loc=cbind(sample_data$km_utm.x, sample_data$km_utm.y),
                    max.edge = c(1,5)*max.edge,
                    cutoff = .5,
                    offset = c(max.edge, -bound.outer))
## Make sure mesh projection matches that of other objects

ggplot() +
  gg(mesh) 
mesh$n



#### SIMULATE BATHY THAT MOVES TO THE RIGHT A BIT

### SIMULATE space1 THAT CHANGES A BIT


library(raster)
library(scales)
library(tidyverse)

logit = function(x) log(x/(1-x))
inv.logit = function(x) exp(x)/(1+exp(x))

dims = 101
mat = matrix(1,ncol=dims,nrow=dims)
n = dims^2

##########################
##### gaus sptial field sim #####
######################
yticks <- xticks <-seq(-3, 3, length=dims)
side <- diff(yticks[1:2])  # side length of squares

######################
### Covariate ######
##################
min_covar=0;max_covar=30
means_cov = 10 + cumsum(seq(.05,1.5,length.out=n_size))
sd_cov = round(seq(4,8,length.out=n_size),1)

####################
### LFD ##############
########################
pop_age_means <- c(1000, 750, 500, 100)
cv=.15
sd_pops = pop_age_means*cv

pop_age_sim = rnorm(4,pop_age_means,sd_pops) ### for st simulation
pop_age_sim;pop_age_means

LFD = c(rnorm(pop_age_sim[1],10,1) , 
        rnorm(pop_age_sim[2],15,2) , 
        rnorm(pop_age_sim[3],25,4), 
        rnorm(pop_age_sim[4],30,8))
hist(LFD,breaks=20)
### discretize in n_size
b = density(LFD)$x
a = seq(min(LFD),max(LFD),length.out=n_size)
dens = density(LFD)$y[sapply(a, function(a, b) {which.min(abs(a-b))}, b)]
barplot(dens,a)


### simulate correlated spatial fields for size
x.k <- book.rspde(coordin, range = range / params[2],
                  sigma = sqrt(params[1]), n = n_size, mesh = mesh,
                  return.attributes = F)
x <- x.k
for (j in 2:n_size) {x[, j] <- rho_size * x[, j - 1] + sqrt(1 - rho_size^2) * x.k[, j]}

scaling_size = scales::rescale(dens,to=c(0,1))

space_effect=space_effect1= effect_covar =  y_Occur = y_count =mu_count=mu_Occur= list()
for(i in 1:n_size){

  covar_profile = round(scales::rescale(4*coordinates(coordin)[,1]-.03*coordinates(coordin)[,1]^2,to=c(min_covar,max_covar)) )+1## same in Y axis
  covar_profile = round(scales::rescale(5*coordinates(coordin)[,1]+coordinates(coordin)[,1]^3,to=c(min_covar,max_covar)) )+1## same in Y axis
  
  plot(coordinates(coordin)[,1],covar_profile)
  
  effect_covar[[i]] =  scales::rescale(dnorm(covar_profile,means_cov[i],sd_cov[i]), 
                                to=c(-1,scaling_size[i]))
  
  space_effect[[i]] = scales::rescale(x[,i],to=c(-1,scaling_size[i]))
  scales::rescale()
  
  
  ### occurrence
  intercept_Occur = dens[i]*5
  nu_Occur = intercept_Occur + effect_covar[[i]] + space_effect[[i]]
  mu_Occur[[i]] = inv.logit(nu_Occur)
  y_Occur[[i]] = rbinom(n, size = 1, prob = mu_Occur[[i]])
  
  ### count
  intercept_count = dens[i]*5
  nu_count = intercept_count + effect_covar[[i]] + space_effect[[i]]
  mu_count[[i]] = exp(nu_count)
  y_count[[i]] = rpois(n, lambda = mu_count[[i]])
  
}

space_effect=as.data.frame(x)
names(space_effect)=paste0("Size_",1:n_size)
names(effect_covar)=paste0("Size_",1:n_size)
names(y_Occur)=paste0("Size_",1:n_size)
names(y_count)=paste0("Size_",1:n_size)
names(mu_Occur)=paste0("Size_",1:n_size)
names(mu_count)=paste0("Size_",1:n_size)

t(lapply(effect_covar,mean))
dens*5
t(lapply(space_effect,mean))



s_effect = as.data.frame(space_effect) %>% 
  mutate(lon = coordinates(coordin)[,1],
         lat = coordinates(coordin)[,2],
         idx=1:n()) %>% 
  pivot_longer(cols= 1:n_size , values_to = "SpatialEffect",names_to = "Size")%>% 
  mutate(size_num = as.numeric(str_sub(Size,6,-1)))

ggplot(s_effect)+geom_point(aes(x=lon,y=lat,color=SpatialEffect))+facet_wrap(~size_num)

c_effect = as.data.frame(effect_covar) %>% 
  mutate(lon = coordinates(coordin)[,1],
         lat = coordinates(coordin)[,2],
         idx=1:n(),
         covar = covar_profile) %>% 
  pivot_longer(cols= 1:n_size , values_to = "CovarEffect",names_to = "Size")%>% 
  mutate(size_num = as.numeric(str_sub(Size,6,-1)))

ggplot(c_effect)+geom_line(aes(x=covar,y=CovarEffect,color=factor(size_num)),linewidth=1)+theme_bw() + 
  #ggtitle("Simulated covariate effects") + 
  xlab("Covariate") + ylab("Effect") +
  theme_bw() + guides(color = guide_legend(title="Size"))
ggsave(filename = "~/bathy_size_sims.png")
ggplot(c_effect)+geom_point(aes(x=covar,y=CovarEffect))+facet_wrap(~Size)

ggplot(c_effect)+geom_tile(aes(x=lon,y=lat,fill=covar)) +theme_bw()+
  xlab("Longitude")+ylab("Latitude")+ #ggtitle("Simulated covariate space")+
  viridis::scale_fill_viridis()+ guides(fill = guide_legend(title="Covariate"))
ggsave(filename = "~/covariate_space.png")


occur = as.data.frame(y_Occur) %>% 
  mutate(lon = coordinates(coordin)[,1],
         lat = coordinates(coordin)[,2],
         idx=1:n()) %>% 
  pivot_longer(cols= 1:n_size , values_to = "Occurrence",names_to = "Size")%>% 
  mutate(size_num = as.numeric(str_sub(Size,6,-1)))

mu_occur = as.data.frame(mu_Occur) %>% 
  mutate(lon = coordinates(coordin)[,1],
         lat = coordinates(coordin)[,2],
         idx=1:n()) %>% 
  pivot_longer(cols= 1:n_size , values_to = "Mean Occur",names_to = "Size")%>% 
  mutate(size_num = as.numeric(str_sub(Size,6,-1)))


count = as.data.frame(y_count) %>% 
  mutate(lon = coordinates(coordin)[,1],
         lat = coordinates(coordin)[,2],
         idx=1:n()) %>% 
  pivot_longer(cols= 1:n_size , values_to = "Count",names_to = "Size") %>% 
  mutate(Count_occur = ifelse(Count>0, 1,0))

mu_count = as.data.frame(mu_count) %>% 
  mutate(lon = coordinates(coordin)[,1],
         lat = coordinates(coordin)[,2],
         idx=1:n()) %>% 
  pivot_longer(cols= 1:n_size , values_to = "Mean Abundance",names_to = "Size") 


SimData = left_join(s_effect,c_effect) %>% left_join(occur) %>% left_join(count) %>% 
  left_join(mu_occur) %>% left_join(mu_count) %>% 
  mutate(size_num = as.numeric(str_sub(Size,6,-1)))



###############################
####### create survey ########
#########################
samp = sample(unique(count$idx),400)
Survey = SimData %>% filter(idx %in% samp)

################
### plots ####
############
SimData %>% filter(size_num==1) %>% ggplot(aes(x=lon,y=lat,color=covar)) + geom_point() 
SimData %>% group_by(Size) %>% ggplot(aes(x=lon,y=lat,color=CovarEffect)) + geom_point() + facet_wrap(~Size)





ggplot(SimData) + geom_boxplot(aes(x=factor(covar), y = Count))+facet_wrap(~Size)
ggplot(SimData) + geom_boxplot(aes(x=factor(covar), y = jitter(Occurrence))) + facet_wrap(~Size)


ggplot(SimData)+ geom_tile(aes(x=lon,y=lat,fill=Count))+facet_wrap(~Size) + 
  geom_point(data=Survey,aes(x=lon,y=lat),color="red")+
  viridis::scale_fill_viridis()
ggplot(SimData)+ geom_tile(aes(x=lon,y=lat,fill=factor(Count_occur)))+facet_wrap(~Size) + geom_point(data=Survey,aes(x=lon,y=lat))
ggplot(SimData)+ geom_tile(aes(x=lon,y=lat,fill=factor(Occurrence)))+facet_wrap(~Size) + geom_point(data=Survey,aes(x=lon,y=lat))

ggplot(SimData)+ geom_tile(aes(x=lon,y=lat,fill=SpatialEffect))+facet_wrap(~Size) + geom_point(data=Survey,aes(x=lon,y=lat),color="red")+
  viridis::scale_fill_viridis()
ggplot(SimData)+ geom_tile(aes(x=lon,y=lat,fill=CovarEffect))+facet_wrap(~Size) + geom_point(data=Survey,aes(x=lon,y=lat),color="red")+
  viridis::scale_fill_viridis()



###################
#### modeling ###
##################
library(inlabru)
library(INLA)
covar_seq = seq(min(SimData$covar),max(SimData$covar),length.out = 16)
mapper_covar <- bru_mapper( INLA::inla.mesh.1d( covar_seq , boundary="free") , indexed = FALSE)


survey_sf <- sf::st_as_sf(Survey, coords = c("lon","lat"))

max.edge = dims*.05
bound.outer = .3
mesh = inla.mesh.2d(loc.domain = cbind(c(0,dims,dims,0,0),c(0,0,dims,dims,0)),
                    #loc=cbind(sample_data$km_utm.x, sample_data$km_utm.y),
                    max.edge = c(1,5)*max.edge,
                    cutoff = .5,
                    offset = c(max.edge, -bound.outer))
## Make sure mesh projection matches that of other objects

ggplot() +
  gg(mesh) +
  #gg(SLDF_UTM, color = "DarkGreen") +
  gg(survey_sf, color = "red")
mesh$n

matern <- inla.spde2.pcmatern(mesh, 
                              prior.range = c(20, 0.1), ## min spatial correlation range, the associated probability that the actual value is below this.
                              prior.sigma = c(0.5, 0.5)) ## max standard deviation, the associated probability that the actual value is greater than this.


cmp_both = Count ~ Intercept(1) + 
  covar(covar, model = "rw2", mapper =mapper_covar, scale.model=T,
        hyper = list(prec = list(prior = "pc.prec", param = c(4, 0.01))),
        group = size_num ,
        ngroup=n_size,
        control.group = list(model="ar1")) +                   
  # size(size_num, model = "rw2", scale.model=T,
  #       hyper = list(prec = list(prior = "pc.prec", param = c(4, 0.01)))) +
  size(Size, model = "factor_contrast") +
  S_size(geometry, model=matern, 
            group = size_num , 
            ngroup = n_size,
            control.group = list(model="ar1"))

fit_both = bru(components = cmp_both,
                 data = survey_sf ,
                 family ="nbinomial",
                 options = list(verbose=T,control.inla = list(int.strategy = "eb")))

plot(fit_both,"size")
plot(fit_both,"covar")

summary(fit_both)

SimData %>% group_by(size_num) %>% summarise(mean=mean(Count)) %>% ggplot(aes(x=size_num,y=mean)) +geom_point()

#### plot mean size in the study area
ggplot(fit_both$summary.random$size)+geom_point(aes(x=ID,y=mean))+  
  geom_errorbar(aes(x=ID,y=mean,ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9)) +
  geom_point(aes(x=1:9,y=rescale(dens,to=c(-.5,1.75))),color="red",size=2) + theme_bw()

summary(fit_both)

fit_both$summary.hyperpar


save(fit_both,file = "~/fit_both.Rdata")
save(SimData,file = "~/SimData.Rdata")
save(survey_sf,file = "~/survey_sf.Rdata")
save(mesh,file = "~/mesh.Rdata")


#####################
##### plots #######
#################

df_space = data.frame(Lat = rep(mesh$loc[,2],n_size),Lon = rep(mesh$loc[,1],n_size),
                mean = fit_both$summary.random$S_size$mean,
                group_size = rep(sort(unique(survey_sf$size_num)),each=mesh$n))

ggplot(df_space) + geom_point( aes(x=Lon, y=Lat, color=mean) ) + facet_wrap(~ group_size) +
  xlim(0,100) + ylim(0,100)

ggplot(SimData)+ geom_tile(aes(x=lon,y=lat,fill=SpatialEffect))+facet_wrap(~Size)+
  viridis::scale_fill_viridis()


df_depth = data.frame(depth = fit_both$summary.random$covar$ID, 
                      estimate =fit_both$summary.random$covar$mean,
                      size = rep(sort(unique(survey_sf$size_num)),each=length(covar_seq)))

ggplot(df_depth) + geom_line(aes(x=depth,y=estimate,color=factor(size)),linewidth=1) + 
  #ggtitle("Estimated covariate effect") + 
  xlab("Covariate") + ylab("Estimate") +
  theme_bw() + guides(color = guide_legend(title="Size"))
ggsave(filename = "~/bathy_size_estimate.png")


#################
### Predict ####
###############
SimData_sf = sf::st_as_sf(SimData, coords = c("lon","lat"))
ea = SimData_sf %>% filter(idx%in%sample(unique(SimData_sf$idx),2000))

predictions <- predict(fit_both, newdata = ea,
                       formula= ~exp(Intercept + covar  + S_size + size),
                       n.samples=2000,
                       include = c("Intercept" , "covar"  , "S_size", "size"))

predictions$lon = st_coordinates(predictions)[,1]
predictions$lat = st_coordinates(predictions)[,2]
ggplot(predictions)+geom_point(aes(x=median,y=Mean.Abundance))+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(~size_num) +
  ggtitle("Predicted vs Simulated mean abundances") + xlab("Predicted") + ylab("Simulated")
ggsave(filename = "~/observed_vs_predicted.png")


################################
#### Size prediction #########
############################
predictions %>% group_by(size_num) %>% summarise(Estimated = sum(median),
                                                 Simulated = sum(Mean.Abundance)) %>%  
  pivot_longer(cols = Estimated:Simulated,names_to = "Origin",values_to = "Abundance") %>% 
  ggplot()+geom_line(aes(x=size_num,y=Abundance,color=Origin),linewidth=1) +#ggtitle("Overall LFD") +
  theme_bw() + xlab("Size")+ylab("Overall abundance")+
  guides(color = guide_legend(title=NULL))
ggsave(filename = "~/Overall_LFD.png")


#################################
#### spatial prediction ########
###############################

# Plot the points using ggplot2
ggplot(predictions) +  # Replace 'X' and 'Y' with your column names
  geom_sf(aes(color=median),size=3) + facet_wrap(~size_num) +#ggtitle("Estimated spatial abundances per size class") +
  viridis::scale_color_viridis()+
  guides(color = guide_legend(title="Estimate"))
ggsave(filename = "~/estimated_maps.png")

ee = as.data.frame(predictions)
ggplot(ee) +  # Replace 'X' and 'Y' with your column names
  geom_point(aes(x=lon,y=lat,color=median),size=4) + facet_wrap(~size_num) +#ggtitle("Estimated spatial abundances per size class") +
  viridis::scale_color_viridis() + xlab("Longitude")+ ylab("Latitude")+
  guides(color = guide_legend(title="Median"))
ggsave(filename = "~/estimated_maps2.png")

ggplot(SimData) +  # Replace 'X' and 'Y' with your column names
  geom_point(aes(x=lon,y=lat,color=Count),size=3) + facet_wrap(~size_num) +#ggtitle("Simulated abundance maps per size class")+
  viridis::scale_color_viridis()+ xlab("Longitude")+ ylab("Latitude")+
  guides(color = guide_legend(title="Abundance"))
ggsave(filename = "~/simulated_maps.png")

ee %>% mutate(Residual = Mean.Abundance-median) %>% ggplot +
  geom_point(aes(x=lon,y=lat,color=Residual),size=4) + facet_wrap(~size_num) +#ggtitle("Spatial residuals per size class")+
  viridis::scale_color_viridis()+ xlab("Longitude")+ ylab("Latitude")+
  guides(color = guide_legend(title="Residual"))
ggsave(filename = "~/residual_maps.png")

##################################
#### covar prediction ##########
##############################
predictions %>% group_by(covar,size_num) %>% summarise(Estimate = mean(median)) %>% 
  ggplot()+geom_point(aes(x=covar,y=size_num,color=Estimate,size=Estimate))+theme_bw() +
  viridis::scale_color_viridis() +#ggtitle("Estimated mean abundance per Size at covariate") +
  theme_bw() + xlab("Covariate")+ ylab("Size")
ggsave(filename = "~/size_covar_mean_map_prediction.png")


predictions %>% group_by(covar,size_num) %>% summarise(Simulated = mean(Mean.Abundance)) %>% 
  ggplot()+geom_point(aes(x=covar,y=size_num,color=Simulated,size=Simulated))+theme_bw()+
  viridis::scale_color_viridis()+#ggtitle("Simulated mean abundance per Size at covariate") +
  theme_bw() + xlab("Covariate")+ ylab("Size")
ggsave(filename = "~/size_covar_mean_map_simulated.png")



predictions %>% group_by(covar,size_num) %>% summarise(Residual = mean(Mean.Abundance-mean)) %>% 
  ggplot()+geom_point(aes(x=covar,y=size_num,color=Residual,size=Residual))+theme_bw()+
  viridis::scale_color_viridis()+#ggtitle("Mean residual per Size at covariate") +
  theme_bw() + xlab("Covariate")+ ylab("Size")
ggsave(filename = "~/size_covar_mean_residual_map.png")



predictions %>% rename(Estimated = "mean",Simulated = "Mean.Abundance") %>% pivot_longer(cols = c("Estimated","Simulated"),names_to = "Type",values_to = "Value") %>%    
  ggplot()+geom_boxplot(aes(x=covar,y=log(Value),color=Type,group=interaction(covar,Type)))+
  facet_wrap(~size_num)+#ggtitle("Estimated abundance at covariate per size class") +
  theme_bw() + xlab("Covariate")+ ylab("Log(Abundance)") + 
  scale_x_continuous(breaks=seq(1, 31, 5))
ggsave(filename = "~/sim_vs_estim_BySize_map.png")


predictions %>% rename(Estimated = "mean",Simulated = "Mean.Abundance") %>% 
mutate(group_covar = cut(covar,breaks=c(1,5,10,15,20,25,31),labels = c("1-5","6-10","11-15","16-20","21-25","26-31"),include.lowest=T) )%>% 
  pivot_longer(cols = c("Estimated","Simulated"),names_to = "Type",values_to = "Value") %>%    
  ggplot()+geom_boxplot(aes(x=size_num,y=log(Value),color=Type,group=interaction(size_num,Type)))+
  facet_wrap(~group_covar)+#ggtitle("Estimated abundance at covariate per size class") +
  theme_bw() + xlab("Size")+ ylab("Log(Abundance)") + 
  scale_x_continuous(breaks=seq(1, 16, 2))
ggsave(filename = "~/sim_vs_estim_CovarGroup_map.png")



predictions %>%   ggplot()+geom_boxplot(aes(x=covar,y=mean,group=covar))+
  facet_wrap(~size_num)+#ggtitle("Estimated abundance at covariate per size class") +
  theme_bw() + xlab("Covariate")+ ylab("Estimated abundance at covariate") + 
  scale_x_continuous(breaks=seq(1, 31, 5))
ggsave(filename = "~/size_covar_prediction.png")

predictions %>%   ggplot()+geom_boxplot(aes(x=covar,y=Mean.Abundance,group=covar))+
  facet_wrap(~size_num)+#ggtitle("Simulated abundance at covariate") +
  theme_bw() + xlab("Covariate")+ ylab("Simulated abundance at covariate") + 
  scale_x_continuous(breaks=seq(1, 31, 5))
ggsave(filename = "~/size_covar_simulated.png")


mean_estimate =predictions %>% group_by(covar,size_num) %>% summarise(mean = mean(mean)) 
predictions %>%   ggplot()+geom_boxplot(aes(x=covar,y=Count,group=covar))+ 
  geom_line(data=mean_estimate,aes(x=covar,y=mean),color="red")+
  facet_wrap(~size_num)+#ggtitle("Simulated abundance and estimated mean abundance at covariate per size class") +
  theme_bw() + xlab("Covariate")+ ylab("Abundance") + 
  scale_x_continuous(breaks=seq(1, 31, 5))
ggsave(filename = "~/size_covar_observedVsMeanFit.png")

mean_estimate =predictions %>% group_by(covar,size_num) %>% summarise(mean = mean(mean)) 
predictions %>%   ggplot()+geom_boxplot(aes(x=covar,y=log(Count),group=covar))+ 
  geom_line(data=mean_estimate,aes(x=covar,y=log(mean)),color="red")+
  facet_wrap(~size_num)+#ggtitle("Simulated abundance and estimated mean abundance at covariate per size class") +
  theme_bw() + xlab("Covariate")+ ylab("Abundance") + 
  scale_x_continuous(breaks=seq(1, 31, 5))
ggsave(filename = "~/size_covar_LOG_observedVsMeanFit.png")



predictions %>% mutate(Residual = Mean.Abundance-mean) %>% 
  ggplot()+geom_boxplot(aes(x=covar,y=Residual,group=covar))+
  facet_wrap(~size_num)+ggtitle("Residuals at covariate per size class") +theme_bw() + xlab("Covariate")+ ylab("Residual") + 
  scale_x_continuous(breaks=seq(1, 31, 5))
ggsave(filename = "~/size_covar_residuals.png")


########################
#### LFDs per region ##
######################

predictions %>%  mutate( Region = case_when(
  st_coordinates(predictions)[, "X"] < 33 & st_coordinates(predictions)[, "Y"] > 66 ~ "Region 1",
  st_coordinates(predictions)[, "X"] >= 33 & st_coordinates(predictions)[, "X"] <= 66 & 
    st_coordinates(predictions)[, "Y"] <= 66 &  st_coordinates(predictions)[, "Y"] >= 33 ~ "Region 2",
  st_coordinates(predictions)[, "X"] >= 66 &  
    st_coordinates(predictions)[, "Y"] <= 33 ~ "Region 3",
  TRUE ~ NA
))   %>%  ggplot()  +  # Replace 'X' and 'Y' with your column names
  geom_sf(aes(color=Region)) + theme_bw() +scale_color_discrete(na.translate=F)
ggsave(filename = "~/spat_regions.png")


predictions %>%  mutate( Region = case_when(
  st_coordinates(predictions)[, "X"] < 33 & st_coordinates(predictions)[, "Y"] > 66 ~ "Region 1",
  st_coordinates(predictions)[, "X"] >= 33 & st_coordinates(predictions)[, "X"] <= 66 & 
    st_coordinates(predictions)[, "Y"] <= 66 &  st_coordinates(predictions)[, "Y"] >= 33 ~ "Region 2",
  st_coordinates(predictions)[, "X"] >= 66 &  
    st_coordinates(predictions)[, "Y"] <= 33 ~ "Region 3",
  TRUE ~ NA
)) %>% filter(!is.na(Region))  %>% group_by(Region) %>% mutate(Total_region_estimated=sum(median),
                                    Total_region_Simulated = sum(Mean.Abundance)) %>% ungroup() %>% 
  group_by(Region,size_num) %>% summarise(Estimated = sum(median/Total_region_estimated),
                                         Simulated = sum(Mean.Abundance/Total_region_Simulated)) %>%  
  pivot_longer(cols = Estimated:Simulated,names_to = "Type",values_to = "Abundance") %>% 
  ggplot()+geom_line(aes(x=size_num,y=Abundance,color=Type),linewidth=1) + 
  facet_wrap(~Region) +theme_bw() + xlab("Size cm") + 
  theme_bw()
ggsave(filename = "~/LFDs_region.png")


########################
#### LFDs per point ##
######################
muestra=sample(unique(predictions$idx),9)

predictions %>% filter(idx%in%muestra) %>% ggplot+geom_sf_text(aes(label = idx), size = 3, color = "black") + theme_bw() +
  xlab("")+ylab("")+ggtitle("LFD at these points in space")
ggsave(filename = "~/spat_points.png")


predictions %>% filter(idx%in%muestra) %>%
  mutate(Estimated = median,
         Simulated = Mean.Abundance) %>% 
  pivot_longer(cols = c(Estimated,Simulated),names_to = "Source",values_to = "Abundance") %>% 
  ggplot()+geom_line(aes(x=size_num,y=Abundance,color=Source),linewidth=1) + facet_wrap(~idx) +theme_bw() + xlab("Size cm") + 
  theme_bw()
ggsave(filename = "~/LFDs_at_point.png")





library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')

dat<- read.csv("Snail Kite Gridded Data_larger.csv", header = T, sep = ",")

#remove IDs w < 200 obs; won't run well
dat<- dat %>% group_by(id) %>% filter(n() > 200) %>% ungroup()
dat.list<- df.to.list(dat=dat)
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval


#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 40000


## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details

dat.res<- behavior_segment(dat = behav.list, ngibbs = ngibbs)
###Takes 73.37 min to run for 40000 iterations for 26 IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)

##Determine MAP for selecting breakpoints
MAP<- apply(dat.res$LML, 1, function(x) getMAP(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, MAP = MAP, brk.cols = 99)  #brk.cols is max matrix cols


## Heatmaps
heatmap(data = behav.list, brkpts = brkpts, dat.res = dat.res, type = "behav")



#########################################
#### Assign Behavioral Time Segments ####
#########################################

dat_out<- map(behav.list, assign.time.seg) %>% map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_LDA_behavior")
write.csv(dat_out, "Snail Kite Gridded Data_larger_behav.csv", row.names = F)



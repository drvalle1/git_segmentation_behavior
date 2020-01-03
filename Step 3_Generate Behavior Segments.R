
library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')

dat<- read.csv("Snail Kite Gridded Data_larger.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()

#if dt within 5 min of 1 hr, round to 1 hr
dat<- round_track_time(dat = dat, int = 3600, tol = 5/60*3600)

dat.list<- df.to.list(dat=dat)
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval
behav.list<- behav.list[sapply(behav.list, nrow) > 2]  #remove IDs w/ fewer than 3 obs


#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 40000


## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details

dat.res<- behavior_segment(dat = behav.list, ngibbs = ngibbs)
###Takes 1 hr to run 40000 iterations for 31 IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, brk.cols = 99)  #brk.cols is max matrix cols


## Heatmaps
plot.heatmap(data = behav.list, nbins = c(6,8), brkpts = brkpts, dat.res = dat.res, type = "behav")



#########################################
#### Assign Behavioral Time Segments ####
#########################################

dat_out<- map(behav.list, assign.time.seg) %>% map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_LDA_behavior")
write.csv(dat_out, "Snail Kite Gridded Data_larger_behav.csv", row.names = F)



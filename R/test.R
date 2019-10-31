require(dplyr)
require(Rcpp)
source("./R/Fourier_est.R")
sourceCpp("./src/Fourier_calc.cpp")

rawdata = read.csv("_data_/CleanEvents.csv",stringsAsFactors = F)
rawdata$Date_Time = as.POSIXct( paste(rawdata$Date,rawdata$Time),format = "%Y/%m/%d %H:%M:%OS")

rawdata = na.omit(rawdata)

mindate_st = min(rawdata$Date_Time[rawdata$Site=="Stockton"]) %>%
  as.character()%>%
  substr(1,10) %>%
  paste("00:00:00") %>%
  as.POSIXct(format = "%Y-%m-%d %H:%M:%OS") %>%
  as.numeric() # first second at the first detection day

rawdata$sec = as.numeric(rawdata$Date_Time)-as.numeric(mindate_st)# all in seconds

# species
spplist = unique(rawdata$Species)

# data of stockton
stockton = lapply(unique(rawdata$Species), function(spp,rawdata){
  rawdata[rawdata$Site=="Stockton" & rawdata$Species==spp,]
},rawdata)

names(stockton) = unique(rawdata$Species)

# time at event
event_time = lapply(stockton,function(kk){sort(kk$sec)})


spplist = names(event_time)[sapply(event_time, length)>3]

event_time_list = lapply(spplist, function(spp,event_time){
  event_time[[spp]]
},event_time)

names(event_time_list) = spplist

period = 24*3600 # period, a day, can be change based on sun


MLE_fourier = lapply(event_time_list,function(ww){
  optim(runif(7),logL,n=3,event_time = ww,method = "BFGS",control = list(maxit = 1000))
}) # MLE using event-time

lapply(MLE_fourier,function(w){w$convergence}) # check convergence

AICs = lapply(MLE_fourier,function(ww){
  2*length(ww$par)+2*ww$value

})

set.seed(42)

MCMC_fourier = lapply(event_time_list,function(ww){
  ActivityPP_sampler(event_time_list$Bear_black,
                     n_sample = 500000,
                     n_burn_in = 50000,
                     thin_by = 250,
                     n=2,n_points = 3000,P=86400,
                     prop_var = .1)
})

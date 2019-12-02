require(dplyr)
require(Rcpp)
source("./R/Fourier_est.R")
sourceCpp("./src/Fourier_calc.cpp")
sourceCpp("./src/STTC.cpp")

rawdata = read.csv("_data_/CleanEvents.csv",stringsAsFactors = F)
rawdata$Date_Time = as.POSIXct( paste(rawdata$Date,rawdata$Date.Time),format = "%Y/%m/%d %H:%M:%OS")

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


spplist = names(event_time)[sapply(event_time, length)>100]

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
  ActivityPP_sampler(ww,
                     n_sample = 50000,
                     n_burn_in = 5000,
                     thin_by = 50,
                     n=2,n_points = 3000,P=86400,
                     prop_var = .1)
})


Time_range = c(0,max(unlist(event_time)))

run_sttc(Time_range,event_time_list$Bear_black,event_time_list$Squirrel,200)

dts = 10^(seq(3,6,0.1))

STTC_dt = sapply(dts,function(dt,Tr,e1,e2){
  run_sttc(Tr,e1,e2,dt)


},Time_range,event_time$Bobcat
,event_time$Fox_red)


plot(log(dts),STTC_dt)

require(ggplot2)


ggplot(data = data.frame(dt=log(dts/3600),STTC = STTC_dt),aes(x=dt,y=STTC))+
  geom_point() +
  geom_line()+
  geom_smooth() +
  geom_hline(yintercept = 0)



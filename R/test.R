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

predic_fourier = lapply(MLE_fourier,function(res){
  lambda_Fourier_series(seq(0,86400,1),3 # number of terms
                        ,res$par[1] # intercept
                        ,res$par[2:4] # amplitude
                        ,86400 # period, a day
                        ,res$par[5:7]) # phase
}) # predict one day's intensity function using the fitted Fourier series

detection_rate = lapply(MLE_fourier,function(res){
  1-exp(-Lambda_Numeric_int(86400,3
                        ,res$par[1]
                        ,res$par[2:4]
                        ,86400
                        ,res$par[5:7],10000))
}) # detection rate per day calculated using integral of detection intensity for a day


plot(predic_fourier$Coyote/(predic_fourier$Fox_red+predic_fourier$Coyote))
# this calculate at a certain time point, given we see a coyote or a fox
#  , what is the probability that is a coyote rather than a fox


test_sample = ActivityPP_sampler(event_time_list$Coyote,n_sample = 100000,n_burn_in = 50000,thin_by = 50,n=2)




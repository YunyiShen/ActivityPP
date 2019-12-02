STTC_ST = function(spike1, spike2, dt, r, time_range, study_area, N_MC=1e6){
	MC = as.data.frame( spsample(study_area,N_MC,"random"))
	MC$t = runif(N_MC,time_range[1],time_range[2])
	
	
	spike1 = as.data.frame(spike1)
	spike2 = as.data.frame(spike2)
	
	N1 = nrow(spike1)
	N2 = nrow(spike2)
	
	entries = c("t","x","y")
	if(prod(entries %in% names(spike1) , entries %in% names(spike2))==0) 
		stop("Names of spike (detection) data should be t, x and y\n")
	
	
	sttc_st = run_sttc(spike1 , N1 , spike2 , N2 , MC , N_MC , dt , r)
	
	return(sttx_st)
}
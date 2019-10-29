logL = function(par,event_time,n=3,P=86400,n_points = 3000){
  A0 = par[1]
  par = par[-1]
  A = par[1:n]
  phi = par[1:n + n]

  -logLCpp(event_time,n,A0,A,P,phi,n_points)
}

ActivityPP_sampler = function(event_time,n_sample = 1000,n_burn_in=1000,thin_by = 1,n=3,P=86400,n_points=3000,prop_var = 0.1){
	require(coda)
	res = mcmc(matrix(NA,nrow = ceiling( n_sample/thin_by),ncol = 2*n+1),thin = thin_by)

	ini = c( runif(n+1),sapply(1:n,function(i,P){
	  runif(1,0 , pi)
	},P) )
	#curr_par = optim(ini,logL,event_time=event_time,n=n,P=P,n_points = n_points,method = "BFGS",control = list(maxit = 1000))$par
	curr_par = ini
	curr_posterior = -logL(curr_par,event_time,n,P,n_points)
	#cat(curr_posterior,"/n")
  prop_par = curr_par
	## NEED TO FIX THE RANGE OF PHASE PARAMETER!!
	for(i in 1:(n_sample+n_burn_in)){
		prop_par[0:n+1] = curr_par[0:n+1] + rnorm(n+1,0,prop_var)
    prop_par[1:n + n + 1] = runif(n,0,pi)
		prop_posterior = -logL(prop_par,event_time,n,P,n_points)
    #cat(prop_posterior,"\n")
		#if(exp(prop_posterior-curr_posterior)>0) cat(exp(prop_posterior-curr_posterior),"\n")
    rn = runif(1)
    #cat(rn,"\n")
    cat("curr:",curr_posterior,"prop:",prop_posterior,"\n")
		if(rn<exp(prop_posterior-curr_posterior)){
			curr_par = prop_par
			curr_posterior = prop_posterior
			cat(i,"accepted\n")
		}




		if(i>n_burn_in & (i-n_burn_in)%%thin_by==0) res[(i-n_burn_in)/thin_by,]=curr_par
	}
	return(res)
}

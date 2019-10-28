Intensity = function(event_time,period,npoints=1001){
  waitings_s =
      c(event_time[1],event_time[-1] - event_time[-length(event_time)])
  npoints = npoints + 1
  Zi = event_time - period * floor(event_time/period) # normalized time to event
  a_hat = mean(floor(waitings_s/period))/(mean(floor(waitings_s/period))+1)# MOM of probability that encounter one in the period

  scaled_Zi = 2*pi*Zi/period
  bw = getBandWidth(scaled_Zi)
  ts = 2*pi*seq(from=0,to=period,length.out = npoints)/period

  den = densityFit(scaled_Zi,ts,bw)

  den = 2*pi*den/period
  interval = period/(npoints-1)
  CDF_hat = interval * sapply(1:(npoints-1),function(i,den){
    sum(den[1:i])
  },den) # numeric integral

  lambda_x = ((1-a_hat)*den[1:npoints-1])/(1-(1-a_hat)*(CDF_hat))

  return(list(t=ts[-npoints],lambda = lambda_x))
}


densityFit <-
  function(x, grid, bw) {
    n <- length(x)
    nxpts <- length(grid)
    dens <- densRad( (x), as.integer(n),(nxpts),
               (grid), (bw))
    return(dens)
  }

getBandWidth <-
  function(A, kmax = 3) {
    estkappa <- numeric(kmax)
    for (k in 1:kmax) {
      trigmom <- TrigMomRad(A, k)
      # CHECK Afun value at end points: opposite signs or zero
      # Upper limit reduced from 500 to 345 in version 0.3.1.9003, 2018-05-03
      if (Afun(0.0001, trigmom, k) * Afun(345, trigmom, k) <= 0)
        estkappa[k] <- uniroot(Afun, c(0.0001,345), trigmom, k)$root
    }
    kappahat <- max(estkappa[1:kmax])
    if(kappahat > 0)  {
      return( ( 3 * length(A) * kappahat^2 * besselI(2 * kappahat, 2) /
                  (4 * sqrt(pi) * besselI(kappahat, 0)^2) )^(2/5) )
    }  else  {  # Would happen if all k gave uniroot errors.
      return(NA)
    }
  }


TrigMomRad <-
  function(x, p) {
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    circmean <- atan2(sinr, cosr)
    sin.p <- mean(sin(p * (x - circmean)))
    cos.p <- mean(cos(p * (x - circmean)))
    sqrt(sin.p^2 + cos.p^2)
  }

Afun <-
  function(kappa, trigmom, k) {
    besselI(kappa, k) / besselI(kappa, 0) - trigmom
  }


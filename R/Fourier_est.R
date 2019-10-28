logL = function(par,event_time,n=3,P=86400,n_points = 3000){
  A0 = par[1]
  par = par[-1]
  A = par[1:n]
  phi = par[1:n + n]
  
  -logLCpp(event_time,n,A0,A,P,phi,n_points)
}


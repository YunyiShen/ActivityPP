#include <Rcpp.h>
using namespace Rcpp;

/// Following intensity function
Rcpp::NumericVector lambda_Fourier_terms(const NumericVector & ts, //nodes evaluated
                            const int & w, // which term
                            const double & A, // amplitude
                            const double & P, // period, known
                            const double & phi){//phase
  Rcpp::NumericVector res;
  res = A * cos((2 * M_PI * w * ts / P) - phi);
  return(res);
}

// [[Rcpp::export]]
Rcpp::NumericVector lambda_Fourier_series(const NumericVector & ts, //nodes evaluated
                                  const int & n, // how many terms
                                  const double & A0, // intercept
                                  const NumericVector & A, // amplitude
                                  const double & P, // period, known
                                  const NumericVector & phi){//phase
  NumericVector res;
  res = 0*ts + A0/2; // first term
  for(int i=0;i<n;i++){
    res = res + lambda_Fourier_terms(ts,i+1,A[i],P,phi[i]);
  }
  return(exp(res)); // make sure large than 0
}


/// take integral at final point, used in likelihood
// [[Rcpp::export]]
double Lambda_Numeric_int(const double & t_last, //nodes evaluated
                                          const int & n, // how many terms
                                          const double & A0, // intercept
                                          const NumericVector & A, // amplitude
                                          const double & P, // period, known
                                          const NumericVector & phi,
                                          const int & n_points){//phase
  double res;
  NumericVector ts;
  ts = seq(0,n_points);
  ts = (P/n_points) * ts;// points for one period
  double interval = P/(n_points-1);
  int n_period = floor(t_last/P); // how may periods at t_last
  double last_int = t_last - n_period * P;
  int n_points_last = floor( n_points * last_int/P); // how may points to take
  NumericVector last_ts = head(ts,n_points_last); // last period, get part of the series

  res = n_period * sum(interval * lambda_Fourier_series(ts,n,A0,A,P,phi)) +
    sum(interval * lambda_Fourier_series(last_ts,n,A0,A,P,phi));

  return(res);
}



// Likelihood
// [[Rcpp::export]]
double logLCpp(const NumericVector & event_time,
            const int & n, // how many terms
            const double & A0, // intercept
            const NumericVector & A, // amplitude
            const double & P, // period, known
            const NumericVector & phi,
            const int & n_points){

  return(sum(log(lambda_Fourier_series(event_time,n,A0,A,P,phi)))
         - Lambda_Numeric_int( event_time[event_time.length()-1],n,A0,A,P,phi,n_points)
  ); // log likelihood function for event_time

}


// MCMC sampler
// [[Rcpp::export]]
NumericMatrix ActivityPP_samplerCPP(const NumericVector & event_time,
								    const int & n_sample,
								    const int & n_burn_in,
								    const int & thin_by,
								    const int & n,
								    const double & P,
								    const int & n_points,
								    const double & prop_var){
	NumericMatrix res(ceil(n_sample/thin_by),2*n+1);

	NumericVector curr(2*n+1);
	curr[Range(0,n)]=runif(n+1);
	curr[Range(n+1,2*n)] = M_PI * runif(n);
	NumericVector prop = curr;

	double curr_post;
	curr_post = logLCpp(event_time,n,curr[0],curr[Range(1,n)],P,curr[Range(n+1,2*n)],n_points);
	//printf("%f\n",((curr_post)));
  double prop_post = 0;
  int r = 0;
  double rn;

	for(int i = 0; i<(n_sample+n_burn_in); ++i){
    prop = curr;
		prop[Range(0,n)] = prop[Range(0,n)]+rnorm(n+1,0,prop_var);
		prop[Range(n+1,2*n)] = runif(n,0,M_PI);
		prop_post = logLCpp(event_time,n,prop[0],prop[Range(1,n)],P,prop[Range(n+1,2*n)],n_points);
		printf("curr:%f , prop:%f \n",(curr_post),prop_post);
    rn=runif(1,0,1)[0];
    //printf("%f\n",(rn));
		if((rn)<exp(prop_post-curr_post)){
		  printf("%d accepted \n",i+1);
			curr = prop;
			curr_post = prop_post;
		}
		//printf("%f\n",curr_post);
		if(((i+1)>n_burn_in) & ( (i+1-n_burn_in)%thin_by==0)){
			r = (i+1-n_burn_in)/thin_by-1;
		  for(int j = 0;j<2*n+1;++j){
		    res(r,j)=curr[j];
		  }

		}

	}

	return(res);


}







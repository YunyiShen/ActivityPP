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


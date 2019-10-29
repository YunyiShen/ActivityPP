#include <Rcpp.h>
using namespace Rcpp;

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

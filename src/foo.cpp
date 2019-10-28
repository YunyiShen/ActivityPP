#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::NumericVector densRad(const Rcpp::NumericVector & x, const int & n,
                   const int & nxpts,
                   const Rcpp::NumericVector & xpts,
                   const double & h,
                   Rcpp::NumericVector res){
  double d, ksum, cons;
  cons = 2 * M_PI * R::bessel_i(h, 0, 2);
  for(int i=0; i < nxpts; i++) {
    ksum = 0;
    for(int j=0; j < n; j++) {
      d = xpts[i] - x[j];
      ksum += pow(exp(cos(d) - 1), h);
    }
    res[i] = ksum / n / cons;
  }
  return(res);
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//



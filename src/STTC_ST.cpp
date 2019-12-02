#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

////
// Code modified from reference:
// Cutts, Catherine S., and Stephen J. Eglen. "Detecting pairwise correlations in spike trains: an objective comparison of methods and application to the study of retinal waves." Journal of Neuroscience 34.43 (2014): 14288-14303.
// https://github.com/CCutts/Detecting_pairwise_correlations_in_spike_trains
////



double dist(double x1,double y1,double x2,double y2){
	return(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)));
}



// function for parallel MC calculating T

struct get_P : public Worker
{
   // source vector
    const RVector<double> spike1_t;
    const RVector<double> spike1_x;
    const RVector<double> spike1_y;

    const int N_spike1;

    const RVector<double> spike2_t;
    const RVector<double> spike2_x;
    const RVector<double> spike2_y;

    const int N_spike2;

	const double dt;
	const double r;

   // accumulated value
   int value;

   // constructors
   get_P(const NumericVector spike1_t,
         const NumericVector spike1_x,
         const NumericVector spike1_y,
         const int N_spike1,
         const NumericVector spike2_t,
         const NumericVector spike2_x,
         const NumericVector spike2_y,
         const int N_spike2,
		 const double dt,
		 const double r) : spike1_t(spike1_t) , spike1_x(spike1_x) , spike1_y(spike1_y) , N_spike1(N_spike1),
                            spike2_t(spike2_t) ,spike2_x(spike2_x) , spike2_y(spike2_y) , N_spike2(N_spike2) , dt(dt) , r(r) ,value(0) {}
   get_P(const get_P& get_P_sample, Split) : spike1_t(get_P_sample.spike1_t) , spike1_x(get_P_sample.spike1_x) , spike1_y(get_P_sample.spike1_y) ,
   							N_spike1(get_P_sample.N_spike1),
                            spike2_t(get_P_sample.spike2_t) ,spike2_x(get_P_sample.spike2_x) , spike2_y(get_P_sample.spike2_y) ,
							N_spike2(get_P_sample.N_spike2) , dt(get_P_sample.dt) , r(get_P_sample.r) , value(0) {}

   // accumulate just the element of the range I've been asked to
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
	  	int flag = 0;
		double spike1_t_temp = spike1_t[i];
		double spike1_x_temp = spike1_x[i];
		double spike1_y_temp = spike1_y[i]; // the MC sample currently work with
        for (std::size_t j = 0; j < N_spike2; ++j) {
			double spike2_t_temp = spike2_t[j];
			double spike2_x_temp = spike2_x[j];
			double spike2_y_temp = spike2_y[j];

			if (spike2_t_temp < spike1_t_temp+dt &&
				spike2_t_temp > spike1_t_temp-dt &&
				dist(spike2_x_temp,spike2_y_temp,spike1_x_temp,spike1_y_temp)<r) {
				flag = 1;
				break; // if it is in a range of some spike, break and flag=1
			}

		 }

	    value += flag;
	  }
   }

   // join my value with that of another get_T
   void join(const get_P& rhs) {
      value += rhs.value;
   }
};



//calculates the tiling coefficient

double run_P(int N1, int N2, double dt, double dr, const List& spike_1, const List& spike_2){
    int j = 0;
    int Nab = 0;
    NumericVector spike_time_1 = spike_1["t"];
    NumericVector spike_x_1 = spike_1["x"];
    NumericVector spike_y_1 = spike_1["y"];

    NumericVector spike_time_2 = spike_2["t"];
    NumericVector spike_x_2 = spike_2["x"];
    NumericVector spike_y_2 = spike_2["y"];

    double r_diff = 0;
    double t_diff = 0;// time and spatial difference between two spikes


    for(int i = 0; i < N1; ++i){
        while(j < N2){
//check every spike in train 1 to see if there's a spike in train 2 within dt  (don't count spike pairs)
// don't need to search all j each iteration

            t_diff = abs(spike_time_1[i] - spike_time_2[j]);
            r_diff = sqrt((spike_x_1[i] - spike_x_2[j]) * (spike_x_1[i] - spike_x_2[j]) + (spike_y_1[i] - spike_y_2[j]) * (spike_y_1[i] - spike_y_2[j]));
            if( t_diff <= dt && r_diff <= dr){
                Nab = Nab+1;
                break;
            }
            else if(spike_time_2[j] > spike_time_1[i]+dt){
                break;// skip all spikes that later than dt of the one we are working with
            }
            else{
                ++j;
            }
        }
    }
    return Nab;
}


// [[Rcpp::export]]
double run_Ppara(int N1, int N2, double dt, double r, const List& spike1, const List& spike2){
    int j = 0;
    int Nab = 0;
    NumericVector spike1_t = spike1["t"];
    NumericVector spike1_x = spike1["x"];
    NumericVector spike1_y = spike1["y"];

    NumericVector spike2_t = spike2["t"];
    NumericVector spike2_x = spike2["x"];
    NumericVector spike2_y = spike2["y"];


	get_P get_p(spike1_t,spike1_x,spike1_y,N1,spike2_t,spike2_x,spike2_y,N2,dt,r);

   // call parallel_reduce to start the work
    parallelReduce(0, N1, get_p);


    return (get_p.value/(double)N1);
}





// [[Rcpp::export]]
double run_sttc(const List& spike1,
				const int& N1,
				const List& spike2,
				const int& N2,
				const List& MC,
				const int& N_MC,
                double dt,
			    double r){

    double TA;
    double TB;
    double PA;
    double PB;
    double index;


	PA = run_Ppara(N1,N2,dt,r,spike1,spike2);
	PB = run_Ppara(N2,N1,dt,r,spike2,spike1);

	TA = run_Ppara(N_MC,N1,dt,r,MC,spike1);
	TB = run_Ppara(N_MC,N2,dt,r,MC,spike2);

    index=0.5*(PA-TB)/(1-TB*PA)+0.5*(PB-TA)/(1-TA*PB);


    return(index);
}

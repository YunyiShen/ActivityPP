#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

////
// Code modified from reference:
// Cutts, Catherine S., and Stephen J. Eglen. "Detecting pairwise correlations in spike trains: an objective comparison of methods and application to the study of retinal waves." Journal of Neuroscience 34.43 (2014): 14288-14303.
// https://github.com/CCutts/Detecting_pairwise_correlations_in_spike_trains
////



//calculates the tiling coefficient
// [[Rcpp::export]]
double run_P(int N1,int N2,double dt,const NumericVector& spike_times_1, const NumericVector& spike_times_2){
	int j = 0;
	int Nab = 0;

	for(int i = 0; i < N1; ++i){
		while(j < N2){
//check every spike in train 1 to see if there's a spike in train 2 within dt  (don't count spike pairs)
// don't need to search all j each iteration
			if(abs(spike_times_1[i] - spike_times_2[j]) <= dt){
				Nab = Nab+1;
				break;
			}
			else if(spike_times_2[j] > spike_times_1[i]){
				break;
			}
			else{
				++j;
			}
		}
	}
	return Nab;
}


// [[Rcpp::export]]
double run_T(int N,double dt,double start, double end, const NumericVector& spike_times){

	//double dt = dtv;
	//double start = startv;
	//double end = endv;
	//int N1 = N1v;
	double time_A;
	int i = 0;
	double diff;

//maximum
	time_A = 2 * (double)N * dt;

// if just one spike in train
	if(N==1){

	  if((spike_times[0] - start) < dt){
	    	time_A = time_A - start+spike_times[0] - dt;
	  }
	  else if((spike_times[0] + dt) > end){
	   	 time_A=time_A - spike_times[0] - dt + end;
	      }

	}

//if more than one spike in train
	else{


			while(i<(N-1)){

				diff = spike_times[i+1] - spike_times[i];

				if(diff < 2*dt){
					//subtract overlap
					time_A = time_A - 2 * dt + diff;

				}

				++i;
			}

			//check if spikes are within dt of the start and/or end, if so just need to subract
			//overlap of first and/or last spike as all within-train overlaps have been accounted for


			if((spike_times[0] - start) < dt){

			  time_A = time_A-start+spike_times[0] - dt;
			}


			if((end - spike_times[N-1]) < dt){

			  time_A = time_A - spike_times[N - 1] - dt + end;
			}
              	}

	return time_A;
}


// [[Rcpp::export]]
double run_sttc(const NumericVector& Time,// first element for starting time and last for end
			  const NumericVector& spike_times_1,// first time-to-event
			  const NumericVector& spike_times_2,
			  double dt){
	int N1 = spike_times_1.length();
	int N2 = spike_times_2.length();
	double TA;
	double TB;
	double PA;
	double PB;
	double T;
	double index;



	if(N1==0 || N2==0){
		return( INT_MIN);
	}
	else{
		T=Time[1]-Time[0];
		TA=run_T(N1,dt,Time[0],Time[1], spike_times_1);
		TA=TA/T;
		TB=run_T(N2,dt,Time[0],Time[1], spike_times_2);
		TB=TB/T;
		PA=run_P(N1,N2,dt, spike_times_1, spike_times_2);
		PA=PA/(double)N1;
		PB=run_P(N2,N1,dt, spike_times_2, spike_times_1);
		PB=PB/(double)N2;
		index=0.5*(PA-TB)/(1-TB*PA)+0.5*(PB-TA)/(1-TA*PB);
	}

	return(index);
}

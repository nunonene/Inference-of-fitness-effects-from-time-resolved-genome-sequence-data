/* 

 emission.cpp

 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 
*/

//own headers...

#include "emission.h"

using namespace std;

// Constructor
Emission::Emission(){
  is_set        = 0; // class is set flag
  mode          = 1; // 1: Multinomial, 2:Dirichlety-multinomial (to be implemented)
  optmode       = 1;
  nSim          = 0; // total number of simulations;
  maxNsim       = 0; // maximum number of simulations to be used in likelihood optimization;
  nTs           = 0; // Total number of sampling instants;
  KH            = 0; // Number of Haplotypes;
  SmplDt        = 0; // Sampling period;
  C             = 100; // Sequencing depth;
  T             = 2000; // Total number of generations;
  N             = 1000000000; // Population size;
  mu            = 0.00001; // Mutation rate;
  DSMod         = 0; // model to be used;
  sTimes        = NULL; // sampling times vector;
  reads         = NULL; // data matrix;
  Type          = 0; // 0:Lin 1:Hyper


}

// Set
void Emission::set(vector<int>& nts){

  sTimes  = new int [nTs];
	

  for (int s=0; s<nTs; s++){
	  sTimes[s]   = nts[s];
  }

  reads  = new unsigned int ** [KH];

   for (int k=0;k<KH;k++){
	   reads[k]  = new unsigned int * [maxNsim];
		   for (int n=0; n<maxNsim; n++){
			   reads[k][n]    = new unsigned int [nTs];
		   }
   }

  is_set=1; // class is set flag
}

// Destructor

Emission::~Emission(){

	if (is_set==1){
		Emission::clear();
	}
}

void Emission::clear(){
  if (is_set==0) abort();


  for (int k=0; k<KH; k++){
    for (int n=0; n<maxNsim; n++){

      delete [] reads[k][n];
    }

    delete reads[k];
  }

  delete [] reads;

  delete [] sTimes;


  is_set   = 0;
}

// Multinomial

double Emission:: multinomial_lnpdf (const size_t K, const double p[], const unsigned int n[])
{
   size_t k;
   unsigned int C = 0;
   double log_pdf = 0.0;
   double norm = 0.0;

   for (k = 0; k < K; k++)
     {
       C += n[k];
     }

   for (k = 0; k < K; k++)
     {
       norm += p[k];
     }

   log_pdf = gsl_sf_lnfact (C);

   for (k = 0; k < K; k++)
     {
       /* Handle case where p[k]==0 and n[k]==0 */


	   if (!(p[k] == 0 && n[k]==0))
         {
           log_pdf += log (p[k] / norm) * n[k] - gsl_sf_lnfact (n[k]);
         }
     }

   return log_pdf;
}

bool value_comparer(std::map<int,double>::value_type &i1, std::map<int,double>::value_type &i2){
  return i1.second < i2.second;
}

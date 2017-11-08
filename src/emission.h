/*
 emission.h
 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 
*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_psi.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_cdf.h"

using namespace std;

class Emission{
public:
  Emission();
  ~Emission();
  void set(vector<int>& stimes);
  double multinomial_lnpdf (const size_t K, const double p[], const unsigned int n[]);
  void clear();
  int mode;
  int optmode;
  int is_set;
  unsigned int *** reads; // vector with sampled data
  unsigned int C; // sequencing depth
  int N; // population size
  double mu; // global mutation rate
  int nSim; // number of simulations or replicates per paramter set
  
  int maxNsim; // maximum number of simulations to be used
  int nTs; // total number of sampling instances
  int SmplDt; // sampling period
  int * sTimes; // sampling generations
  int T; // Total number of generations
  int KH; // number of haplotypes
  int DSMod;
  int Type;
};

bool value_comparer(std::map<int,double>::value_type &i1, std::map<int,double>::value_type &i2);

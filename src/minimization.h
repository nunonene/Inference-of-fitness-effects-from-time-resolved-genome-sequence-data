/*

  minimization.h
  
  Created on: Nov 2017
  
  Author: Nuno Rocha Nene
  E-mail:  nunonene@gmail.com
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#include <list>


// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_linalg.h"
#include <gsl/gsl_sort_vector.h>

// Own headers...

#include "emission.h"
#include "PropMod.h"

using namespace std;

// Optimization struct

struct fpar{ // structure of parameters used for optimization...

  PropMod * myProp;
  int to_opt;
  double s_i;
  double s_f;
  int sgrid;
  gsl_matrix * range;
  double scurr;
};


// The  numerical minimization routine to get selection parameters only with MC method (used for likelihood optimization under a deterministic or a stochastic model)...

double find_local_optimum_MC_sel (gsl_rng *rgen, gsl_vector * other, gsl_matrix * range,
		  void * params,
		  double (*obj_fn)( const gsl_vector * x, void * p),
		  double prec,
		  int& steps,
		  int verbose);


// The  numerical minimization routine to get selection and threshold parameters with MC method (used for likelihood optimization under delay-deterministic models)...

double find_local_optimum_MC (gsl_rng *rgen, gsl_vector * other, gsl_matrix * range,
		  void * params,
		  double (*obj_fn)( const gsl_vector * x, void * p),
		  double prec,
		  int& steps,
		  int verbose);

// The  numerical minimization routine to get parameters by a hierarchical MC method...(only used for likelihood optimization under a K beta delay-deterministic model)...

double find_local_optimum_MChier (gsl_rng *rgen, gsl_vector * other, gsl_matrix * range,
		  void * params,
		  double (*obj_fn)( const gsl_vector * x, void * p),
		  double prec,
		  int& steps,
		  int verbose);


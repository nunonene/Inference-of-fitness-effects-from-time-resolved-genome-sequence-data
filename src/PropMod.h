/*
 
 PropMod.h

 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 
*/

#ifndef SRC_PROPMOD_H_
#define SRC_PROPMOD_H_

#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <map>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// GSL headers...

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multimin.h"


using namespace std;



class PropMod{

public:
  PropMod();
  ~PropMod();
  Emission * myEmit;
  int is_set;
  int Rep;
  
  int mode;
  double s;
  double * beta;

  double ** betag;

  double * nThresh;
  gsl_matrix * MutProp;
  gsl_vector * FSPL;
  gsl_vector * FTR;

  gsl_rng *rgen;

  int Nstch;



  double ** mean_ML;

  double *** mean_ML_STOCH;

  vector<int> is_identity;

  double total_llh;

  double * llhsim;

  double do_Fwd(int s);

  void set(Emission * emit);

  void setRep(Emission * emit);

  void setMutMat(gsl_matrix *& MProp);

  void setFSPL(gsl_vector *& SPL);

  void setFTR(gsl_vector *& FR);

  void setRep(int rp);

  void clear();

  int Det();
  
  int DetT1();

  int TRESH();

  int TRESHglobal();

  int StochABC();

  double get_total_llh();

  double get_llh();

};



#endif /* SRC_PROPMOD_H_ */

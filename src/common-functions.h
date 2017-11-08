/*

 common-functions.h

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
#include <algorithm>
#include <ctime>


// gsl headers

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"



// Own headers

class Emission;

using namespace std;

int get_dims(const char * data_fn,int& maxNsim, int& Khaplo,int& T); // Check if dimensions of data collected from the file are cosnsistent from those provided from command line
void get_data(const char * data_fn, Emission * myEmit, vector<int>& selectlines); // import data

void get_mutgraph(gsl_matrix * MutProp,gsl_vector * FSPL,int Type,int KH); // Get weighted mutation graph for the linear fitness chain or the fitness hypercube

void get_HyperGraphMutMat(gsl_matrix * MutProp,int n); // import mutation network weights

void get_HyperGraphSPL(gsl_vector * FSPL,int n,int Type); // import distance between haplotypes


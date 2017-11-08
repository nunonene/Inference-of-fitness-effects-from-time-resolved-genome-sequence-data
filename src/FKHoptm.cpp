/*
 
 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 


 Code used in paper 

 "A delay-deterministic model for inferring fitness effects from time-resolved genome sequence data", Nuno R. Nene, Alastair S. Dunham,  Christopher J. R. Illingworth. (submitted)


 for inference of fitness landscapes from simulated Wright-Fisher evolutionary time-series.


 Additional code for inference from experimental time-series is available elsewhere (see paper for details).

 FKHoptm.cpp

  Output:

  -- Inferred fitness values s and additional delay-deterministic model threshold (\beta) parameters. Population size N is not inferred.


  Required Arguments:

  --data      [str]     Observed time-frequencies file (each line corresponds to a replicate) when read depth was C (see below)

  --out       [str]     Output file

  --C         [int]     Sequencing read depth (simulated)

  --N         [int]     Population size (provided not estimated)

  --mu        [double]  Mutation rate (here, known)

  --Nsim      [int]     Number of replicates (trajectories)

  --maxNsim   [int]     Number of simulations to be used (if absent = Nsim)
 
  --T         [int]     Length of evolutionary trajectories, i.e. number of generations.

  --nTs       [int]     Number of sampling instances to be used

  --Khaplo    [int]     Number of haplotypes in the simulations (each haplo is in a different file)

  --si        [double]  Initial selection amplitude (default minimum value si=0)

  --sf        [double]  Final   selection amplitude (default minimum value sf=2)

  --sgrid     [int]     Number of initial points for parameter optimization

  --DSModel   [int]     0: Deterministic model (without knowledge of N);
			1: Delay-deterministic model with mutation Heaviside function, one extra parameter \betai per node (without knowledge of N);
		        2: Stochastic model (with knowledge of N);
		        3: Delay-deterministic model with mutation Heaviside function, one per node with global extra parameter \beta (without knowledge of N);
		        4: Delay-deterministic model with mutation Heaviside function, one per node with global extra parameter \beta (without knowledge of N) and 			           information on survival probability \pi

  --mode      [int]     Emission model: 1 Multinomial

  --optmode   [int]     [1,3]     Optimization routine; 3 is to be used only with DSModel = 1;

  --AllSim    [0,1]     0: Use all replicates for estimates 1: Perform optimization on a replicate by replicate basis

  --Type      [0,1]     0: Linear Fitness Chain; 1: Fitness Hypercube


 Notes:

 1. In the paper the initial conditions were known; at t=0 all the frequency was concentrated at the fitness minimum
 
 2. --data should be in the form of Nsim replicates by T columns; there should be one file per haplotype and each should be named *_Haplo*_smpl.txt.

 3. Files with mutational weights for the hypercube should be HyperGraphMutMat_nh*.txt, where * represents the number of haplotypes. Any other weights    distribution should be in this format. The sum of out mutation weights should sum to 1. See examples in HyperGraphMutFit folder

 Example (see data in Examples folder):

 (Stochastic model for linear fitness chain)

 ./FKHoptm --data SimKhaploDynamics_Khaplo6_N1000000000_Nsim100_C100_Init1_sigma0.5_mu0.00001_T2000_TypeLin --out  SimKhaploDynamics_Khaplo6_N1000000000_Nsim100_C100_Init1_sigma0.5_mu0.00001_T2000_TypeLin --N 1000000000 --C 100 --mu 0.00001 --Nsim 100 --maxNsim 100 --T 2000 --nTs 200 --Khaplo 6 --si 0.01 --sf 1 --sgrid 100 --DSModel 2 --mode 1 --optmode 1 --Type 0;

 
*/

// headers

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

//own headers...

#include "minimization.h"
#include "common-functions.h"


using namespace std;



// Command line options structure

struct cmdl_opts{
  int mode;     // Emission model mode
  int optmode;  // Optimization mode
  double seed;  // Seed used in optimization routine
  double s;     // Fitness value (parameter to be optimized)
  double s_i;   // Starting initial condition for fitness s on the grid of s's
  double s_f;   // Final initial condition for fitness s on the grid of s's
  int sgrid;    // Size of grid for fitness s optimization routine
  int Npop;     // Population size;
  int C;        // Sequencing read depth;
  int T;        // Duration of experiment in generations;
  int nTs;      // Number of collected sampling points, here assumed to be through sequencing
  int Nsim;     // Number of simulations starting at the initial value;
  int maxNsim;  // Maximum number of simulations to be used;
  const char *out; // Output filename;
  const char *data_fn; // Data filename;
  double mu;    // Mutation rate (here, known);
  int Khaplo;   // K haplotypes;
  int Type;     // 0 (Linear Fitness chain) or 1 (Fitness Hypercube);
  int DSMod;    // 0: Deterministic model (without knowledge of N);
		// 1: Delay-deterministic model with mutation Heaviside function, one extra parameter \betai per node (without knowledge of N);
		// 2: Stochastic model (with knowledge of N);
		// 3: Delay-deterministic model with mutation Heaviside function, one per node with global extra parameter \beta (without knowledge of N);
		// 4: Delay-deterministic model with mutation Heaviside function, one per node with global extra parameter \beta (without knowledge of N) and 			      information on survival probability \pi
  int nvar;     // Number of parameters to optimize
};


// Maximum likelihood parameter structure

struct MLparams{
	double * llhsim_ML; // maximum Likelihood for each Replicate
	double ** beta_ML; // ML paramters threshold
	double ** n_ML; // ML paramters non-linearity
	double * s_ML;// Fitness paramter
};


// Own Functions


void get_opts( int argc, const char ** argv, cmdl_opts& opts); // get command line options

void test_opts(cmdl_opts& opts); // test dimensions

void default_opts(cmdl_opts& opts); // get default options

void ClearParams (MLparams& MLps,cmdl_opts& opts);

double Q( const gsl_vector * x, void * p); // function for gsl optimization

double find_parameters(gsl_rng * rgen,PropMod * myProp, cmdl_opts& opts); // execute optimization

void init_parameter(gsl_vector* var, gsl_matrix* range, int to_opt, PropMod * myProp, cmdl_opts& opts); // set initial parameters



// *** MAIN***

int main (int argc, const char * argv[]){

    cmdl_opts opts; //  parameters struture

    MLparams MLps; //   auxiliary structure for parameters being optimized

    // Get options

    default_opts(opts);


    get_opts( argc, argv, opts);

    //srand(opts.seed);

    vector<int> sTimes; // vector of sampling instances for each haplotype


    // Get dimensions

    int KH;

    KH=opts.Khaplo;

    const char *file;

    file=opts.data_fn;

    int maxT=get_dims( file, opts.Nsim, KH,opts.T); // check if dims read from files are consistent with parameters read from the command line
    

    // Announce:

    printf("\n ============================================================================= \n");
    printf("\n Fitness Landscape Inference From Partially Observed Wright-Fisher Time-Series \n");
    printf("\n ============================================================================= \n");

    printf("\n--- Hidden model: ");

    if(opts.DSMod == 0){

    	printf("Deterministic \n");

    } else if(opts.DSMod == 1){

        printf("K beta \n");
   
    } else if(opts.DSMod == 2){

        printf("Stochastic \n");
   
    } else if(opts.DSMod == 3){

        printf("1 beta \n");
   
    } else if(opts.DSMod == 4){

        printf("1 beta pi \n");
   
    }

    printf("\n--- Emission model: ");

    if(opts.mode == 1){

    	printf("Multinomial \n");

    }

    if(opts.Type == 1){

    	printf("\n--- Type: fitness hypercube \n");

    }else{

	printf("\n--- Type: linear fitness chain \n");
    }

    printf("\n--- Number of haplotypes: %i\n", KH);

    printf("\n--- Number of trajectories used in the inference of parameters: %i\n", opts.maxNsim);

    printf("\n--- Length of time-series: %i\n ", opts.T);
   

    // Set sampling times

    if(maxT<opts.nTs)opts.nTs=maxT;

    int dt=floor((double)(maxT)/(double)(opts.nTs));


    if(dt==0){ 

	cout<<"Number of sampling generations has to be smaller than the maximum number of generations used: maxT="<<maxT;

	exit(1);
    }

    int count=0;

    int nSmplTs=0;

    
   // Get sampling instants
 	
    for(int i=0;i<maxT;i++){
	if(i==0){
    		sTimes.push_back(i);
    		// cout<<"Sample:"<<i<<" "<<sTimes[nSmplTs]<< "\n";
		 nSmplTs++;
    	}else if(count==dt){
    		sTimes.push_back(i);
    		// cout<<"Sample:"<<i<<" "<<sTimes[nSmplTs]<< "\n";
    		 nSmplTs++;
    		count=0;
    	}
    	count++;
    }

    // Set Emission class model

    Emission myEmit;

    myEmit.T               = maxT; // Total number of generations to be used

    myEmit.mode            = opts.mode; // Emission mode 1: Multinomial 2: Dirichlet-multinomial (to be implemented)

    myEmit.optmode         = opts.optmode; // Optimization mode: 0-Simplex 1:Simple MC

    myEmit.C               = opts.C; // Sequencing depth

    myEmit.DSMod           = opts.DSMod; // Propagation model

    myEmit.nTs             = (int)sTimes.size(); // Number of sampling instances

    myEmit.nSim            = opts.Nsim; // Total number of trajectories/replicates

    myEmit.maxNsim         = opts.maxNsim; // Maximum number of trajectories/replicates used in paramter inference

    myEmit.KH              = KH; // Number of haplotypes

    myEmit.Type            = opts.Type; // 0: Linear chain 1: Hypercube

    myEmit.mu              = opts.mu; // Mutation rate

    myEmit.N               = opts.Npop; // Population size



    myEmit.set(sTimes); // set sampling/sequencing times/generations


    // Get data and input into Emit class

    vector<int> selectlines;

    get_data(file, &myEmit,selectlines);


    // Get weighted mutation graph

    gsl_matrix * MutProp = gsl_matrix_alloc (KH,KH); // mutation weights matrix to be used in propagation

    gsl_vector * FSPL = gsl_vector_alloc (KH); // vector of distances from all haplotypes to fitness maximum is simulated landscape

    gsl_matrix_set_zero(MutProp);

    get_mutgraph(MutProp,FSPL,opts.Type,KH); // import (hypercube) or create (linear) weighted mutation graph


   // Set gsl random generator


    PropMod myProp; //Propagation class

    int s= opts.seed;

    gsl_rng_env_setup();

    gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);

    gsl_rng_set (rgen, s);


    if(opts.DSMod==2){ // Necessary for stochastic model

    	myProp.rgen = gsl_rng_alloc (gsl_rng_taus);

    	gsl_rng_set (myProp.rgen,s);
    }


    // Initialize ML parameter structure


	MLps.llhsim_ML   = new double [opts.maxNsim]; // likelihood    		
	MLps.s_ML   = new double [opts.maxNsim]; // fitness values

    	if(opts.DSMod==3 || opts.DSMod==4){

    		MLps.beta_ML   = new double * [opts.maxNsim];

    		for (int t=0; t<opts.maxNsim; t++){

    		    MLps.beta_ML[t]   = new double [1];

    		 }
    	}


    	if(opts.DSMod==1){
    		MLps.beta_ML   = new double * [opts.maxNsim];
    		 
		for (int t=0; t<opts.maxNsim; t++){
    		     MLps.beta_ML[t]   = new double [KH];
    		 }
    	}


	
    // Find maximum-likelihood estimates of all parameters in each model for each replicate ...

    	for (int t=0; t<opts.maxNsim; t++){

    		if(t==0){
    			myProp.set( &myEmit);

    			myProp.setMutMat(MutProp);

    			gsl_matrix_free(MutProp);

    			myProp.setFSPL(FSPL);

    			gsl_vector_free(FSPL);

    			myProp.setRep(t);
    		}else{

    			myProp.setRep(t);
    		}

		MLps.llhsim_ML[t] = find_parameters(rgen, &myProp, opts); // optimize
    				
    		MLps.s_ML[t]=myProp.s;


    		for (int k=0; k<KH; k++){
    			if(opts.DSMod==1){
    				MLps.beta_ML[t][k]=myProp.beta[k];
    			}
    			
    			if(opts.DSMod==3 || opts.DSMod==4){
    				MLps.beta_ML[t][0]=myProp.beta[0];
    			}
    		}
    	}
    

    gsl_rng_free (myProp.rgen);


    gsl_rng_free (rgen); // free gsl_random generator


   
    // Export ML parameters and respective likelihoods

    char buffparams[1024];


    if(opts.DSMod==0){

	sprintf(buffparams,"%s.paramsMLDet_nTs%i.txt", opts.out,opts.nTs);

    } else if (opts.DSMod==1){

        if(opts.optmode==3){
		sprintf(buffparams,"%s.paramsML1betaKbeta_nTs%i.txt",  opts.out,opts.nTs);			
	}else{
		sprintf(buffparams,"%s.paramsMLKbeta_nTs%i.txt", opts.out,opts.nTs);				
	}

    } else if (opts.DSMod==2){
      
	sprintf(buffparams,"%s.paramsMLStoch_nTs%i.txt", opts.out,opts.nTs);
	
    } else if (opts.DSMod==3){
	sprintf(buffparams,"%s.paramsML1beta_nTs%i.txt", opts.out,opts.nTs);


    } else if (opts.DSMod==4){

	sprintf(buffparams,"%s.paramsML1betaPi_nTs%i.txt", opts.out,opts.nTs);

    }
     
    //
    
    FILE * params_estimates = fopen(buffparams,"w"); // open ML parameters file
    	
    for (int t=0; t<opts.maxNsim; t++){
    	  
         fprintf(params_estimates, "%d %f %f ",selectlines[t],MLps.s_ML[t],MLps.llhsim_ML[t]);
    		   
        if(opts.DSMod==1){
    		      
	   for (int i=1;i<opts.nvar;i++){
							  
		fprintf(params_estimates, "%f ",MLps.beta_ML[t][i-1]);
           }
				 
	 }else if(opts.DSMod==3 || opts.DSMod==4){

	    for (int i=1;i<opts.nvar;i++){
								 
		fprintf(params_estimates, "%f ",MLps.beta_ML[t][0]);
	    }	
	}

	fprintf(params_estimates, "\n");
     }
      

    // Clear structures

       ClearParams(MLps,opts);

       selectlines.clear();

    // Close file

       fclose(params_estimates);

	cout << "\n Results are summarized in:\n\n"<< buffparams << "\n\n";

    return (0);
}

// *** MAIN END ***


// Default options 

void default_opts(cmdl_opts& opts){
    opts.data_fn         = NULL            ; // Data file;
    opts.out             = NULL            ; // Output data file;
    opts.Npop            = 1000000000      ; // Population size;
    opts.mu              = 0.00001         ; // Mutation rate;
    opts.C               = 100             ; // Depth
    opts.sgrid           = 5               ; // Test starting seed values of fitness s at 5 points from sigma_i to sigma_f;
    opts.s               = -1.0            ; // Fitness (not to be changed);
    opts.s_i             = 0.001           ; // seed value for starting PropMod amplitude tested in the likelihood optimization function;
    opts.s_f             = 1               ; // seed value for starting PropMod amplitude tested in the likelihood optimization function;
    opts.sgrid           = 100              ; // number of starting points for PropMod parameter optimization;
    opts.mode            = 1               ; // default option is the Binomial emission model;(other options will be added in future versions)
    opts.optmode         = 1               ; // default option 1, other option is 3 (used only with K beta delay-deterministic model)
    opts.seed            = (int) time(NULL);
    opts.nTs             = 1               ; // Sampling period;
    opts.T               = 2000            ; // number of total generations; 
    opts.Nsim            = 100             ; // Maximum number of replicates to be used
    opts.maxNsim         = opts.Nsim       ;  // Maximum number of simulations to be used;
    opts.Type            = 0               ; // 0:Linear chain; 1: Hypercube;
    opts.DSMod           = 0               ; // 0: Deterministic model (without knowledge of N);
			                     // 1: Delay-deterministic model with mutation Heaviside function, one extra parameter \betai per node (without knowledge of N);
		                             // 2: Stochastic model (with knowledge of N);
		                             // 3: Delay-deterministic model with mutation Heaviside function, one per node with global extra parameter \beta (without knowledge of N);
		                             // 4: Delay-deterministic model with mutation Heaviside function, one per node with global extra parameter \beta (without knowledge of N) and 			      information on survival probability \pi
}

// Get command line arguments...

void get_opts( int argc, const char ** argv, cmdl_opts& opts){

  default_opts(opts);

  int opt_idx = 1;

  int Tp=0;
  int NN=0;
  int CC=0;
  int MM=0;
  int SS=0;
  int TT=0;
  int KK=0;
  int DD=0;

  string opt_switch;  

  while ( opt_idx < argc && (argv[opt_idx][0] == '-')){
    
    opt_switch = argv[opt_idx];
    opt_idx++;

    if (opt_idx==argc) break;
    if ( argv[opt_idx][0] == '-') continue;
    if ( opt_switch.compare("--data") == 0){
          opts.data_fn = argv[opt_idx];
    }else if ( opt_switch.compare("--out") == 0){
          opts.out = argv[opt_idx];
    }else if ( opt_switch.compare("--N") == 0){
          opts.Npop = atoi(argv[opt_idx]);
	  NN=1;
    }else if ( opt_switch.compare("--C") == 0){
          opts.C = atoi(argv[opt_idx]);
	  CC=1;
    }else if ( opt_switch.compare("--mu") == 0){
          opts.mu = atof(argv[opt_idx]);
	  MM=1;
    }else if ( opt_switch.compare("--Nsim") == 0){
          opts.Nsim = atoi(argv[opt_idx]);
          SS=1;
    }else if ( opt_switch.compare("--maxNsim") == 0){
          opts.maxNsim = atoi(argv[opt_idx]);
    }else if ( opt_switch.compare("--T") == 0){
          opts.T = atoi(argv[opt_idx]);
          TT=1;
    }else if ( opt_switch.compare("--nTs") == 0){
          opts.nTs = atoi(argv[opt_idx]);
    }else if ( opt_switch.compare("--Khaplo") == 0){
          opts.Khaplo = atoi(argv[opt_idx]);
	  KK=1;
    }else if ( opt_switch.compare("--si") == 0){
          opts.s_i = atof(argv[opt_idx]);
    }else if ( opt_switch.compare("--sf") == 0){
          opts.s_f = atof(argv[opt_idx]);
    }else if ( opt_switch.compare("--sgrid") == 0){
  	  opts.sgrid = atoi(argv[opt_idx]);
    }else if ( opt_switch.compare("--DSModel") == 0){
          opts.DSMod = atoi(argv[opt_idx]);
          DD=1;
    }else if ( opt_switch.compare("--mode") == 0){
          opts.mode = atoi(argv[opt_idx]);
    }else if ( opt_switch.compare("--optmode") == 0){
          opts.optmode = atoi(argv[opt_idx]);
    }else if ( opt_switch.compare("--Type") == 0){
          opts.Type = atoi(argv[opt_idx]);
          Tp=1;
    }else if ( opt_switch.compare("--seed") == 0){
        opts.seed = atoi(argv[opt_idx]);
    }else {

      cout << "Wrong command line arguments\n"<<endl;
      exit(1);

    }
    opt_switch.clear();
    opt_idx++;
  }

  if(NN==0) {
     cout << "Please provide population size: --N [int]\n"<<endl;
  }

  if(CC==0) {
     cout << "Please provide sequencing read depth: --C [int]\n"<<endl;
  }

  if(MM==0) {
     cout << "Please provide mutation rate: --mu [double]\n"<<endl;
  }

  if(SS==0) {
     cout << "Please provide number of replicates in file: --Nsim [int]\n"<<endl;
  }

  if(TT==0) {
     cout << "Please provide evolutionary time-series length: --T [int]\n"<<endl;
  }

  if(KK==0) {
     cout << "Please provide number of haplotypes: --Khaplo [int]\n"<<endl;
  }

  if(DD==0) {
     cout << "Please provide model: --DSModel [0,1,2,3,4]\n"<<endl;
  }

  if(Tp==0) {
     cout << "Please provide Type of mutation network: --Type [0/1]\n"<<endl;
  }

  test_opts(opts); // test for restriction in combinations of parameters
}

// Test some options provided from command line ...

void test_opts(cmdl_opts& opts){
  
  if ( opts.data_fn == NULL){
      cout<<"ERROR: provide data file prefix\n";
      exit(1);
  }

  if ( opts.out == NULL){
        cout<<"ERROR: provide output filename\n";
        exit(1);
  }


  if(opts.optmode==3 && opts.DSMod!=1){ // test if model and optimization mode are compatible

	cout<<"ATTENTION:--optmode 3 has to be used only with --DSModel 1";

	exit(1);
  }	
}


// Find parameters ...

double find_parameters(gsl_rng *rgen,PropMod * myProp, cmdl_opts& opts){

	int to_opt;
	int nvar=0;
	double llh=0;

	if(opts.DSMod==0){
		nvar=1;
		to_opt=0;
	}else if(opts.DSMod==1){
		nvar=opts.Khaplo+1;
		to_opt=1;
	}else if(opts.DSMod==2){
		nvar=1;
		to_opt=2;
	}else if(opts.DSMod==3){
		nvar=2;
		to_opt=3;
	}else if(opts.DSMod==4){
		nvar=2;
		to_opt=4;
	}

	opts.nvar=nvar; // store number of parameters to optimize


	gsl_vector * var;
	gsl_matrix * range;

	gsl_matrix *varglobal;

	gsl_vector *llhsglobal;

    if(nvar>0){

    	varglobal   = gsl_matrix_alloc(opts.maxNsim,nvar);

	gsl_matrix_set_zero(varglobal);

	llhsglobal   = gsl_vector_alloc(opts.maxNsim);

	gsl_vector_set_zero(llhsglobal);


    	range = gsl_matrix_alloc(nvar,2);
	var   = gsl_vector_alloc(nvar);

    	// Set initial values
        init_parameter(var, range,  to_opt, myProp, opts);



        fpar myfpar;


	myfpar.myProp      = myProp;
        myfpar.to_opt      = to_opt;
        myfpar.s_i         = opts.s_i;
        myfpar.s_f         = opts.s_f;
        myfpar.sgrid       = opts.sgrid;
        myfpar.range=gsl_matrix_alloc(nvar,2);

        gsl_matrix_memcpy ( myfpar.range, range);


        void * param = static_cast<void*>(&myfpar);
    
        // Get the ML estimates and ML value

        int steps = 0;
      
        if(myProp->myEmit->optmode==1){

        	if(myProp->myEmit->DSMod==0 || myProp->myEmit->DSMod==2 || myProp->myEmit->DSMod==3 || myProp->myEmit->DSMod==4){

        		llh = - find_local_optimum_MC_sel(rgen, var, range,  param, &Q, 1.0e-3, steps, 1);
        	}else{

        		llh = - find_local_optimum_MC(rgen, var, range,  param, &Q, 1.0e-3, steps, 1);
        	}

        }else if(myProp->myEmit->optmode==3){

        	llh = - find_local_optimum_MChier(rgen, var, range,  param, &Q, 1.0e-3, steps, 1); // only used for K beta model by pre-computing a 1 beta model
        }


      

        //set the ML values into the objects

        if(to_opt == 0 || to_opt == 2){
                myProp->s = var->data[0];
        }else if(to_opt == 1){
            	myProp->s = var->data[0];
            	for (int i=1; i<nvar; i++){
            	myProp->beta[i-1] = var->data[i];
            	}
        }else if(to_opt == 3 || to_opt == 4){
            	myProp->s = var->data[0];
            	myProp->beta[0] = var->data[1];
        }

       gsl_vector_free(var);
        
       gsl_matrix_free(range);
          
       gsl_matrix_free(myfpar.range);

    }else{
    	cout<< "ERROR in find_D_parameters(): number of paramters 0"<<"\n";
    }

    return(llh);
}

// Initialize parameters ...

void init_parameter(gsl_vector* var, gsl_matrix* range, int to_opt, PropMod * myProp, cmdl_opts& opts){ // set initial values and parameter range

	int nvar = (int) range->size1;

	 if(to_opt == 0 || to_opt == 2){
		var->data[0]=opts.s_i;
		range->data[0]=0.001;
		range->data[1]=2.0;

	}else if(to_opt == 1){
		var->data[0]=opts.s_i;
		range->data[0]=0.0;
		range->data[1]=1.0;

		for (int i=1; i<nvar; i++){
			var->data[i]=1.0;
			range->data[i*2]=0.0;
			range->data[i*2+1]=1.0;
		}

	}else if(to_opt == 3 || to_opt == 4){
		var->data[0]=opts.s_i;
		range->data[0]=0.0;
		range->data[1]=2.0;

		var->data[1]=1.0;
		range->data[2]=0.0;
		range->data[3]=1.0;
	}
}

// Clear ML parameter structure ...

void ClearParams (MLparams& MLps,cmdl_opts& opts){
		
delete [] MLps.llhsim_ML;		
delete [] MLps.s_ML;

		
if(opts.DSMod==1){

	for (int t=0; t<opts.maxNsim; t++){
		delete [] MLps.beta_ML[t];			
	}
			
delete [] MLps.beta_ML;
		
}
		
if(opts.DSMod==3 || opts.DSMod==4 ){		
	for (int t=0; t<opts.maxNsim; t++){				
		delete [] MLps.beta_ML[t];
	}
	delete [] MLps.beta_ML;		
}

}

// Function to be used with gsl optimization ...

double Q( const gsl_vector * x, void * p){ 

  fpar * myfpar = static_cast<fpar*> (p);
  int nvar;

	  if(myfpar->to_opt==0 || myfpar->to_opt==2){
		 nvar=1;
	  } else if(myfpar->to_opt==3 || myfpar->to_opt==4){
		 nvar=2;
	  } else if(myfpar->to_opt==1){
		 nvar=1+(myfpar->myProp->myEmit->KH);
	  } 


  // Check if params are within limits

	double lower;
	double upper;



  if(myfpar->myProp->myEmit->optmode==0){

	  if(!(myfpar->to_opt==0 || myfpar->to_opt==2)){

		  for (int i=0;i<nvar;i++){

			  int a= myfpar->to_opt==6 ? i+1:i;

			  lower=gsl_matrix_get(myfpar->range,a,0);
			  upper=gsl_matrix_get(myfpar->range,a,1);

			  if(x->data[i] < lower){

				  return(1.0e20);
			
			  }else if(x->data[i] > upper){
				  return(1.0e20);
			  }
		  }
	  }
  }


   if(myfpar->to_opt == 0 || myfpar->to_opt == 2){
	   myfpar->myProp->s = x->data[0];
   }else if(myfpar->to_opt == 1){
	   myfpar->myProp->s = x->data[0];
	   for (int i=1; i<nvar; i++){
		   myfpar->myProp->beta[i-1] = x->data[i];
	   }
   }else if(myfpar->to_opt == 3 || myfpar->to_opt == 4){
	   myfpar->myProp->s = x->data[0];
	   myfpar->myProp->beta[0] = x->data[1];

   }

  // Do Fwd to get LLH
    
    double llh=0.0;

    llh = (myfpar->myProp)->get_llh();
  

  return(-llh);
}


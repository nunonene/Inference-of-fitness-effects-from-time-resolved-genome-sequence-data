
/*

 PropMod.cpp
 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 
*/

//own headers...

#include "emission.h"
#include "PropMod.h"
#include "common-functions.h"

// GSL headers

#include <gsl/gsl_matrix.h>


#define PI 3.1415926

// Constructor

PropMod::PropMod(){

  is_set          = 0; // class is set flag
  myEmit          = NULL;
  Rep             = 0;
  s               = -1.0;
  beta            = NULL; // frequency threshold for delay deterministic model
  betag           = NULL;
  nThresh         = NULL;
  mean_ML         = NULL; // predicted mean
  mean_ML_STOCH   = NULL; // predicted mean
  total_llh       = 0;
  myEmit          = NULL;
  mode            = 1;
  llhsim          = NULL;
  MutProp         = NULL;
  FSPL            = NULL;
  Nstch           = 1;
  FTR             = NULL;
  rgen            = NULL;


}

// Destructor

PropMod::~PropMod(){

  PropMod::clear();

}

// Set PropMod class

void PropMod::set(Emission * emit){

	myEmit    = emit;
	mode      = myEmit->mode;

	Nstch=100; // number of stochastic simulations for the stochastic model

	if(myEmit->DSMod==2){

		mean_ML_STOCH=new double ** [Nstch];

		for (int t=0;t<Nstch;t++){

				mean_ML_STOCH[t]=new double * [myEmit->KH];

				for (int i=0;i<myEmit->KH;i++){

					mean_ML_STOCH[t][i]= new double [myEmit->T];
				}
		}

	}else{
		mean_ML=new double * [myEmit->KH];

		for (int t=0;t<myEmit->KH;t++){
			mean_ML[t]=new double [myEmit->T];
		}

	        llhsim=new double[myEmit->nSim];

		
		if(myEmit->DSMod ==1){

			beta=new double[myEmit->KH];

		}else if(myEmit->DSMod ==3 || myEmit->DSMod ==4){
			
			beta=new double[1];

		}
	}

	is_set=1;
}

// Set mutation matrix

void PropMod::setMutMat(gsl_matrix *& MProp){

	 if (MutProp == NULL){

		 MutProp = gsl_matrix_alloc (myEmit->KH,myEmit->KH);

		 gsl_matrix_memcpy (MutProp, MProp);
	 }

	is_set=1;
}

// Set distance vector for each haplotype

void PropMod::setFSPL(gsl_vector *& SPL){

	 if (FSPL == NULL){

		 FSPL = gsl_vector_alloc (myEmit->KH);

		 gsl_vector_memcpy (FSPL, SPL);
	 }

	is_set=1;
}

void PropMod::setFTR(gsl_vector *& FR){



	 if (FTR == NULL){

		 FTR = gsl_vector_alloc (myEmit->KH);
		 //gsl_matrix_set_zero(MutProp);

		 gsl_vector_memcpy (FTR, FR);
	 }



	is_set=1;
}


// Set current repl if optimization is on a repl by repl basis

void PropMod::setRep(int rp){

	Rep       = rp;

	is_set=1;
}

// Clear class vectors and matrices

void PropMod::clear(){

	if (is_set==1){

		if(beta!=NULL){
			delete [] beta;
		}



		if(nThresh!=NULL){
			delete [] nThresh;
		}



		delete [] llhsim;

		if(myEmit->DSMod==2){

				
			for (int t=0;t<Nstch;t++){
		
				for (int i=0;i<myEmit->KH;i++){
							
					delete [] mean_ML_STOCH[t][i];	
				}
					delete [] mean_ML_STOCH[t];	
			}
				
			delete [] mean_ML_STOCH;

		}else{

			for (int i=0;i<myEmit->KH;i++){

				delete [] mean_ML[i];
			}

			delete [] mean_ML;
		}

	gsl_matrix_free(MutProp);

	}

  is_set   = 0;

}

// Calculate likelihood function value

double PropMod::get_total_llh(){

  int nl;
  total_llh   = 0.0;

  // All simulations in parallel if structure is available:




#ifdef _OPENMP
  int nt = min(myEmit->maxNsim, omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for (nl=0; nl<myEmit->maxNsim; nl++){

    double llh = PropMod::do_Fwd(nl);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      total_llh += llh;
    }

  }


  return(total_llh);
}

double PropMod::get_llh(){

    double llh = PropMod::do_Fwd(Rep);
  return(llh);
}




// Forward step 

double PropMod::do_Fwd(int rp){

	// Prepare fwd...

	double llh=0.0;
	double p=0.0;

	double sumx=0.0;

	unsigned int * nc=new unsigned int[myEmit->KH];

	double * x=new double[myEmit->KH];


	const size_t C=myEmit->C;

	int totalC;
	
	int minc;

	int posminc;

	// Forward Pass



	if(myEmit->DSMod == 0){
	  
		PropMod::Det();
	
	}else if(myEmit->DSMod == 1){
	 
		PropMod::TRESH();

	}else if(myEmit->DSMod == 2){
	  
		PropMod::StochABC();
	
	}else if(myEmit->DSMod == 3){
		  
		PropMod::DetT1();
	
	}else if(myEmit->DSMod == 4){
		
		PropMod::TRESHglobal();
	
	}



	if(myEmit->DSMod == 2){

		

		double pT=0.0;
					
		for (int l=0; l< (myEmit->nTs); l++){

			p=0.0;
			totalC=0;
			minc=2*(int)C;
			posminc=0;

			for (int k=0;k<myEmit->KH;k++){

				nc[k] = myEmit->reads[k][rp][l];

				if(nc[k]<minc) {
					minc=nc[k];
					posminc=k;
				}

				totalC+=nc[k];
			}
						
			if(totalC!=C){
				nc[posminc]+=abs(totalC-(int)C);
			}

			totalC=0;

			for (int k=0;k<myEmit->KH;k++){
				totalC+=nc[k];
			}

			if(totalC!=C){

				cout<<"ERROR-1 in PropMod::do_Fwd\n";
				exit(1);
			}

			sumx=0.0;
			pT=0.0;
							
			for (int numsim=0;numsim<Nstch;numsim++){

				for (int k=0;k<myEmit->KH;k++){

					x[k]=mean_ML_STOCH[numsim][k][myEmit->sTimes[l]];

				}
							
				p = gsl_ran_multinomial_lnpdf (myEmit->KH,x,nc);

				if(!isinf(p)){

					pT=pT+exp(p);
				};
			 }


			pT=pT/(double)Nstch;
			p=log(pT);

			if(p>0){
				cout << "Warning: LLH positive, possible error with probabilities";
			}


			if(isinf(p)){
				p=-1e3;
			}

			if(p!=p){ // if pmf is NaN
				p=-1e-3;
			}

			llh=llh+p;

		}

	}else{


		for (int l=0; l< (myEmit->nTs); l++){

			p=0.0;
			totalC=0;
			minc=0;
			posminc=0;

			for (int k=0;k<myEmit->KH;k++){
				nc[k] = myEmit->reads[k][rp][l];

				if(nc[k]<minc) {
					minc=nc[k];
					posminc=k;
				}

				totalC+=nc[k];

				x[k] = mean_ML[k][myEmit->sTimes[l]];



			}
		
			if(totalC!=C){
				nc[posminc]+=abs(totalC-(int)C);

			}

			totalC=0;

			for (int k=0;k<myEmit->KH;k++){

				totalC+=nc[k];

			}

			if(totalC!=C){

				cout<<"ERROR-1 in PropMod::do_Fwd\n";

				exit(1);

			}

			p = myEmit->multinomial_lnpdf (myEmit->KH,x,nc);

			if(p>0){
				cout<< "Warning: LLH positive, possible error with probabilities";
			}

			if(isinf(p)){
				p=-1e3;
			}

			if(p!=p){// if pmf is NaN
				p=-1e3;
			}

			llh=llh+p;
		}
	}

  delete [] x;
  delete [] nc;
  return(llh);
}


// Deterministic model ...

int PropMod::Det(){

	gsl_matrix * mean=NULL;
	gsl_vector * mean_p=NULL;

	// Set Fitness

	gsl_vector * Fitness=NULL;

	if (Fitness == NULL){

		Fitness = gsl_vector_alloc (myEmit->KH);


		for (int k=0; k<myEmit->KH; k++){

			if(myEmit->Type){ // For the fitness hypercube

				gsl_vector_set(Fitness,k,(1+(gsl_vector_get(FSPL,k)+1)*s));

			}else{ // For the liner fitness chain

				gsl_vector_set(Fitness,k,(1+(k+1)*s)); // s is the current estimate of fitness parameter

			}
		}
	}


	if (mean == NULL){
					
		mean = gsl_matrix_alloc (myEmit->KH,myEmit->T);

		gsl_matrix_set_zero(mean);
	}


	if (mean_p == NULL){

		mean_p = gsl_vector_alloc (myEmit->KH);

		gsl_vector_set_zero(mean_p);
	}


        for (int k=0; k<myEmit->KH; k++){

		if(k==0){
		    	gsl_matrix_set(mean,k,0,1); // In the paper the initial conditions were known; at t=0 all the frequency was concentrated at the fitness minimum
		    	mean_ML[k][0]=1;
		}else{
		    	gsl_matrix_set(mean,k,0,0);
		    	mean_ML[k][0]=0;
		}
        }

	double p;

	vector<double> pp;



	gsl_vector * mem = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre_mean = gsl_vector_alloc(myEmit->KH);




	for(int j=1;j<myEmit->T;j++){

		for(int k=0;k<myEmit->KH;k++){
			p=gsl_matrix_get(mean,k,j-1);
			gsl_vector_set(pre,k,p*(double)(myEmit->mu));
			gsl_vector_set(pre_mean,k,p);
		}
		
		gsl_blas_dgemv(CblasNoTrans, 1, MutProp, pre, 0.0, mem);


		gsl_blas_daxpy (1.0, pre_mean,mem);

		double norm=0.0;

		norm=gsl_blas_dasum (mem);

		gsl_vector_scale(mem,1.0/norm);

		norm=0.0;

		gsl_vector_set_zero(mean_p);

		gsl_vector_memcpy(mean_p,mem);



		gsl_vector_mul(mean_p,Fitness); // Apply fitness

		norm=gsl_blas_dasum (mean_p);

		gsl_vector_scale(mean_p,1.0/norm);

		
		for(int k=0;k<myEmit->KH;k++){
			p=gsl_vector_get(mean_p,k);
			gsl_matrix_set(mean,k,j,p);
			mean_ML[k][j]=p; // Create final mean for LLH calculation
		}

	}
	gsl_vector_free(Fitness);
	gsl_vector_free(mem);
	gsl_vector_free(pre);
	gsl_vector_free(pre_mean);

	gsl_vector_free(mean_p);

	gsl_matrix_free(mean);



	//done
	return (0);


}





// K beta delay-deterministic model

int PropMod::TRESH(){



	gsl_matrix * mean=NULL;
	gsl_vector * mean_p=NULL;

	// Set Fitness
	gsl_vector * Fitness=NULL;

	if (Fitness == NULL){
		Fitness = gsl_vector_alloc (myEmit->KH);
		for (int k=0; k<myEmit->KH; k++){
			if(myEmit->Type){
				gsl_vector_set(Fitness,k,(1+(gsl_vector_get(FSPL,k)+1)*s));
			}else{
				gsl_vector_set(Fitness,k,(1+(k+1)*s)); // s is the current estimate of fitness parameter
			}
		}
	}

	if (mean == NULL){
			mean = gsl_matrix_alloc (myEmit->KH,myEmit->T);
		   	gsl_matrix_set_zero(mean);
	}
	if (mean_p == NULL){
		mean_p = gsl_vector_alloc (myEmit->KH);
		gsl_vector_set_zero(mean_p);
	}


   	for (int k=0; k<myEmit->KH; k++){
		    if(k==0){
		    	gsl_matrix_set(mean,k,0,1);

		    	mean_ML[k][0]=1; // In the paper the initial conditions were known; at t=0 all the frequency was concentrated at the fitness minimum

		    }else{
		    	gsl_matrix_set(mean,k,0,0);
		    	mean_ML[k][0]=0;
		    }
    }

	double p;

	vector<double> pp;



	gsl_vector * mem = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre_mean = gsl_vector_alloc(myEmit->KH);




	for(int j=1;j<myEmit->T;j++){



		for(int k=0;k<myEmit->KH;k++){

			p=gsl_matrix_get(mean,k,j-1);

			
			if(p >= beta[k]){
				gsl_vector_set(pre,k,p*(myEmit->mu));
			}else{
				gsl_vector_set(pre,k,0);
			}


			gsl_vector_set(pre_mean,k,p);

		}
		


		gsl_blas_dgemv(CblasNoTrans, 1, MutProp, pre, 0.0, mem);


		gsl_blas_daxpy (1.0, pre_mean,mem);

		

		double norm=0.0;

		norm=gsl_blas_dasum (mem);

		gsl_vector_scale(mem,1.0/norm);



		norm=0.0;

		gsl_vector_set_zero(mean_p);

		gsl_vector_memcpy(mean_p,mem);



		gsl_vector_mul(mean_p,Fitness); // Apply fitness

		norm=gsl_blas_dasum (mean_p);

		gsl_vector_scale(mean_p,1.0/norm);

		for(int k=0;k<myEmit->KH;k++){
			p=gsl_vector_get(mean_p,k);
			gsl_matrix_set(mean,k,j,p);



			mean_ML[k][j]=p; // Create final mean for LLH calculation
		}

	}
	gsl_vector_free(Fitness);
	gsl_vector_free(mem);
	gsl_vector_free(pre);
	gsl_vector_free(pre_mean);

	gsl_vector_free(mean_p);

	gsl_matrix_free(mean);


	//Done
	return (0);

}

// Stochastic model

int PropMod::StochABC(){



	gsl_matrix * mean=NULL;
	gsl_vector * mean_p=NULL;

	// Set Fitness
	
	double * Fitness=new double[myEmit->KH];

	
	for (int k=0; k<myEmit->KH; k++){

		if(myEmit->Type){
			Fitness[k]=(1+(gsl_vector_get(FSPL,k)+1)*s);
		}else{
			Fitness[k]=(1+(k+1)*s);
		}
	}

	if (mean == NULL){
		mean = gsl_matrix_alloc (myEmit->KH,myEmit->T);
		gsl_matrix_set_zero(mean);
	}
	if (mean_p == NULL){
		mean_p = gsl_vector_alloc (myEmit->KH);
		gsl_vector_set_zero(mean_p);
	}



	for(int numsim=0;numsim<Nstch;numsim++){

		for (int k=0; k<myEmit->KH; k++){

			if(k==0){
		    	gsl_matrix_set(mean,k,0,1);

		    	mean_ML_STOCH[numsim][k][0]=1; // In the paper the initial conditions were known; at t=0 all the frequency was concentrated at the fitness minimum

		    }else{
		    	gsl_matrix_set(mean,k,0,0);
		    	mean_ML_STOCH[numsim][k][0]=0;
		    }
		}
	}

	double p;

	double m;

	double *pp;

	 pp = new double[myEmit->KH];

	unsigned int *n;


	n = new unsigned int[myEmit->KH];

	gsl_vector * mem = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre_mean = gsl_vector_alloc(myEmit->KH);

	double norm=0.0;

	for(int numsim=0;numsim<Nstch;numsim++){

		for(int j=1;j<myEmit->T;j++){

			norm=0.0;

			gsl_vector_set_zero(mem);
			gsl_vector_set_zero(pre);
			gsl_vector_set_zero(pre_mean);


			for(int k=0;k<myEmit->KH;k++){
				p=gsl_matrix_get(mean,k,j-1);

				if(p!=0){

					m =  gsl_ran_poisson (rgen,p*(double)(myEmit->mu)*(double)(myEmit->N));
				}else{

					m=0;
				}

				gsl_vector_set(pre,k,((double)m/(double)(myEmit->N)));

				gsl_vector_set(pre_mean,k,p);

				n[k]=0;

			}
		
			gsl_blas_dgemv(CblasNoTrans, 1, MutProp, pre, 0.0, mem);


			gsl_blas_daxpy (1.0, pre_mean,mem);

			norm=gsl_blas_dasum (mem);

			gsl_vector_scale(mem,1.0/norm);

			norm=0.0;

			gsl_vector_set_zero(mean_p);

		
			for(int k=0;k<myEmit->KH;k++){
			    	p=gsl_vector_get(mem,k);

			    	pp[k]=Fitness[k]*p;

			    	norm+=pp[k];

			    	n[k]=0;
			}


			for(int k=0;k<myEmit->KH;k++){
			    	pp[k]=(double)pp[k]/(double)norm;
			}


			gsl_ran_multinomial (rgen, myEmit->KH,myEmit->N, pp, n); // propagate



			for(int k=0;k<myEmit->KH;k++){
				p=(double)n[k]/(double)(myEmit->N);

				gsl_vector_set(mean_p,k,p);



			}
			norm=gsl_blas_dasum (mean_p);

			gsl_vector_scale(mean_p,1.0/norm);

			for(int k=0;k<myEmit->KH;k++){

				p=gsl_vector_get(mean_p,k);
				gsl_matrix_set(mean,k,j,p);

				if(p!=p)exit(1);

				mean_ML_STOCH[numsim][k][j]=p;

			}

		}
	}


	delete [] Fitness;
	delete [] pp;
	delete [] n;

	gsl_vector_free(mem);
	gsl_vector_free(pre);
	gsl_vector_free(pre_mean);

	gsl_vector_free(mean_p);

	gsl_matrix_free(mean);



	//done
	return (0);

}


// 1 beta delay-deterministic model

int PropMod::DetT1(){



	gsl_matrix * mean=NULL;
	gsl_vector * mean_p=NULL;

	// Set Fitness
	gsl_vector * Fitness=NULL;




	if (Fitness == NULL){


		Fitness = gsl_vector_alloc (myEmit->KH);


		for (int k=0; k<myEmit->KH; k++){
	   		if(myEmit->Type){

				gsl_vector_set(Fitness,k,(1+(gsl_vector_get(FSPL,k)+1)*s));

			}else{
				gsl_vector_set(Fitness,k,(1+(k+1)*s)); // s is the current estimate of fitness parameter

			}
		}
	}


	if (mean == NULL){
			mean = gsl_matrix_alloc (myEmit->KH,myEmit->T);
		   	gsl_matrix_set_zero(mean);
	}
	if (mean_p == NULL){
		mean_p = gsl_vector_alloc (myEmit->KH);
		gsl_vector_set_zero(mean_p);
	}


    for (int k=0; k<myEmit->KH; k++){
		    if(k==0){
		    	gsl_matrix_set(mean,k,0,1); // In the paper the initial conditions were known; at t=0 all the frequency was concentrated at the fitness minimum

		    	mean_ML[k][0]=1;

		    }else{
		    	gsl_matrix_set(mean,k,0,0);
		    	mean_ML[k][0]=0;
		    }
    }

	double p;

	vector<double> pp;



	gsl_vector * mem = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre_mean = gsl_vector_alloc(myEmit->KH);




	for(int j=1;j<myEmit->T;j++){


		for(int k=0;k<myEmit->KH;k++){

			p=gsl_matrix_get(mean,k,j-1);


			if(p >= beta[0]){
				gsl_vector_set(pre,k,p*(myEmit->mu));


			}else{
				gsl_vector_set(pre,k,0);
			}


			gsl_vector_set(pre_mean,k,p);

		}
		
		gsl_blas_dgemv(CblasNoTrans, 1, MutProp, pre, 0.0, mem);


		gsl_blas_daxpy (1.0, pre_mean,mem);

		double norm=0.0;

		norm=gsl_blas_dasum (mem);

		gsl_vector_scale(mem,1.0/norm);



		norm=0.0;

		gsl_vector_set_zero(mean_p);

		gsl_vector_memcpy(mean_p,mem);



		gsl_vector_mul(mean_p,Fitness); // Apply fitness

		norm=gsl_blas_dasum (mean_p);

		gsl_vector_scale(mean_p,1.0/norm);

		
		for(int k=0;k<myEmit->KH;k++){
			p=gsl_vector_get(mean_p,k);
			gsl_matrix_set(mean,k,j,p);
			mean_ML[k][j]=p; // Create final mean for LLH calculation
		}

	}
	gsl_vector_free(Fitness);
	gsl_vector_free(mem);
	gsl_vector_free(pre);
	gsl_vector_free(pre_mean);

	gsl_vector_free(mean_p);

	gsl_matrix_free(mean);



	// Done
	return (0);


}



// 1beta pi delay-deterministic model

int PropMod::TRESHglobal(){



	gsl_matrix * mean=NULL;
	gsl_vector * mean_p=NULL;

	// Set Fitness
	gsl_vector * Fitness=NULL;




	if (Fitness == NULL){


		Fitness = gsl_vector_alloc (myEmit->KH);


		for (int k=0; k<myEmit->KH; k++){
	   		if(myEmit->Type){

				gsl_vector_set(Fitness,k,(1+(gsl_vector_get(FSPL,k)+1)*s));

			}else{
				gsl_vector_set(Fitness,k,(1+(k+1)*s)); // s is the current estimate of fitness parameter

			}
		}
	}


	if (mean == NULL){
			mean = gsl_matrix_alloc (myEmit->KH,myEmit->T);
		   	gsl_matrix_set_zero(mean);
	}
	if (mean_p == NULL){
		mean_p = gsl_vector_alloc (myEmit->KH);
		gsl_vector_set_zero(mean_p);
	}


    for (int k=0; k<myEmit->KH; k++){
		    if(k==0){
		    	gsl_matrix_set(mean,k,0,1);

		    	mean_ML[k][0]=1; // In the paper the initial conditions were known; at t=0 all the frequency was concentrated at the fitness minimum

		    }else{
		    	gsl_matrix_set(mean,k,0,0);
		    	mean_ML[k][0]=0;
		    }
    }

	double p;

	vector<double> pp;



	gsl_vector * mem = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre = gsl_vector_alloc(myEmit->KH);
	gsl_vector * pre_mean = gsl_vector_alloc(myEmit->KH);




	for(int j=1;j<myEmit->T;j++){


		for(int k=0;k<myEmit->KH;k++){

			p=gsl_matrix_get(mean,k,j-1);




			if(p >= ((double)beta[0]/(double)(1-exp(-2*(gsl_vector_get(Fitness,k)-1))))){
				gsl_vector_set(pre,k,p*(myEmit->mu));

			}else{
				gsl_vector_set(pre,k,0);
			}

			gsl_vector_set(pre_mean,k,p);

		}

		gsl_blas_dgemv(CblasNoTrans, 1, MutProp, pre, 0.0, mem);


		gsl_blas_daxpy (1.0, pre_mean,mem);

		double norm=0.0;

		norm=gsl_blas_dasum (mem);

		gsl_vector_scale(mem,1.0/norm);

		norm=0.0;

		gsl_vector_set_zero(mean_p);

		gsl_vector_memcpy(mean_p,mem);



		gsl_vector_mul(mean_p,Fitness); // Apply fitness

		norm=gsl_blas_dasum (mean_p);

		gsl_vector_scale(mean_p,1.0/norm);

		for(int k=0;k<myEmit->KH;k++){
			p=gsl_vector_get(mean_p,k);
			gsl_matrix_set(mean,k,j,p);



			mean_ML[k][j]=p; // Create final mean for LLH calculation
		}

	}
	gsl_vector_free(Fitness);
	gsl_vector_free(mem);
	gsl_vector_free(pre);
	gsl_vector_free(pre_mean);

	gsl_vector_free(mean_p);

	gsl_matrix_free(mean);



	// Done
	return (0);


}


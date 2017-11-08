/* 

 minimization.cpp

 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 
*/


#include "minimization.h"

// The  numerical minimization routine to get selection parameters only with MC method (used for likelihood optimization under a deterministic or a stochastic model)...

double find_local_optimum_MC_sel (gsl_rng *rgen, gsl_vector * other, gsl_matrix * range,
		  void * params,
		  double (*obj_fn)( const gsl_vector * x, void * p),
		  double prec,
		  int& steps,
		  int verbose){

	int iter=0,max_iter =10000;
	
        int min_iter=1000; 

	int status=1,ct=0;

	double f1=1.0e20,f2=0;

	int direction;

	double dxmax=0.5;

	double changepar=dxmax;

	double a_rate=0.0;

	int acceptpar=0;

	int trypar=0;

	double move;

	int nvar = (int) other->size;
	gsl_vector *  xcurr=gsl_vector_alloc(nvar);//allocate here
	gsl_vector *  xlastbest=gsl_vector_alloc(nvar);//allocate here
	gsl_vector *  xbest=gsl_vector_alloc(nvar);//allocate here

	// Sigma grid of starting values

	fpar * myfpar = static_cast<fpar*> (params);

	double si = myfpar->s_i;
	double sf = myfpar->s_f;
	int sgrid = myfpar->sgrid;

	double upper=0.0;

	double lower=0.0;


	double fbest=1.0e20;

	double ftest=0.0;

	// Try on a grid of points 

	for(int tt=0;tt<sgrid;tt++){ // since we are optimizing for one variable only a grid of starting values is doable

		for(int i=0;i<nvar;i++){

			lower=si; //gsl_matrix_get(range,i,0);

			upper=sf; //gsl_matrix_get(range,i,1);

			gsl_vector_set(xcurr,i,lower+(((double)tt/(double)sgrid)*(upper-lower)));

		}

		f2=obj_fn(xcurr,params); // compute likelihood for current parameters

		if(f2<f1) {
			gsl_vector_memcpy(xlastbest,xcurr);
			f1=f2;
		}
	}


	gsl_vector_memcpy(xbest,xlastbest);

	fbest=f1; // store best value so far

	steps=0;

	f2=0;

	// Regular MC ...

	do{

		iter++;
			
		steps++;

		if (iter%100==0&&iter>0) {
			a_rate=(acceptpar+0.)/(trypar+0.);
			changepar=changepar*(0.95+a_rate);
			acceptpar=0;
			trypar=0;
		}

		for (int i=0; i<nvar; i++){

			move=(gsl_rng_uniform(rgen)*changepar)-(changepar/2);

			direction=gsl_rng_get(rgen);

			direction=(direction>0) ? 1:-1;
				
			move=direction*gsl_rng_uniform(rgen)*dxmax;
	
           		gsl_vector_set( xcurr, i,gsl_vector_get(xlastbest,i)+move);

           		lower=gsl_matrix_get(range,i,0);

           		upper=gsl_matrix_get(range,i,1);

           		if(gsl_vector_get(xcurr,i)<lower) gsl_vector_set(xcurr, i,lower+fabs(gsl_vector_get(xcurr,i)-lower));


           		if(gsl_vector_get(xcurr,i)>upper) gsl_vector_set(xcurr, i,upper-fabs(gsl_vector_get(xcurr,i)-upper));	
		}

		trypar++;
         
	        f2=obj_fn(xcurr,params); // compute likelihood for current parameters

	        if(f2<f1){
	            f1=f2;
	            gsl_vector_memcpy(xlastbest,xcurr);
	         }

		// Convergence test

	         if (iter > min_iter && f1 != f2){
	                if (f1-f2 > 0.0 && f1-f2 < prec){
	                    ct++;
	                    if (ct==min_iter) status = 0;
			    f1=f2;
			    gsl_vector_memcpy(xlastbest,xcurr);
	                }
	                else{
	                    ct=0;
	                }

	         } else if(iter > min_iter && f1 == f2){
	            	ct++;

	                if (ct==min_iter){
                           status = 0;	
                        }

	                f1=f2;
			gsl_vector_memcpy(xlastbest,xcurr);
	         }
	}
	while (status &&  iter < max_iter);

	ftest =f1;

	if(ftest<fbest){ // test if new value is better

		gsl_vector_memcpy(xbest,xlastbest);
		fbest=ftest;
	}

	// Store (update) best values ...

	for(int i=0;i<nvar;i++){
		other->data[i]=gsl_vector_get(xbest,i);
	}

	 // Cleanup...

	gsl_vector_free(xcurr);
	gsl_vector_free(xlastbest);
	gsl_vector_free(xbest); 



	return(fbest);

}

// The  numerical minimization routine to get selection and threshold parameters with MC method (used for likelihood optimization under delay-deterministic models)...


double find_local_optimum_MC (gsl_rng *rgen, gsl_vector * other, gsl_matrix * range,
		  void * params,
		  double (*obj_fn)( const gsl_vector * x, void * p),
		  double prec,
		  int& steps,
		  int verbose){



	int iter,max_iter =10000;

	int num_iter_pre=100000;
	int min_iter=1000;
	int status,ct=0;
	double f1=0, f2=1.0e20;

	int direction;



	double move;

	int nvar = (int) other->size;
	gsl_vector *  xcurr=gsl_vector_alloc(nvar);
	gsl_vector *  xcurr_max=gsl_vector_alloc(nvar);
	gsl_vector *  xlastbest=gsl_vector_alloc(nvar);

	// gsl_vector *  xbest = gsl_vector_alloc(nvar);


	//fpar * myfpar = static_cast<fpar*> (params);


	double upper=0.0;

	double lower=0.0;

	//double ftest=0.0;

	double dxmax=0.5;

	double changepar=dxmax;

	double a_rate=0.0;

	int acceptpar=0;

	int trypar=0;

        status=1;

        iter=0;

	steps=0;

	ct=0;

	

	// Preliminary search for a optimum

	for(int tt=0;tt<num_iter_pre;tt++){	
		for(int i=0;i<nvar;i++){
			lower=gsl_matrix_get(range,i,0);
			upper=gsl_matrix_get(range,i,1);
			gsl_vector_set(xcurr,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
		}

		if(tt==0) {
			gsl_vector_memcpy(xlastbest,xcurr);
			f1=obj_fn(xcurr,params); // compute likelihood for current parameters
		}else{
			f2=obj_fn(xcurr,params);
			
		}

		if(f2<f1){
			f1=f2;
			gsl_vector_memcpy(xlastbest,xcurr);
		}		
	}

		
		
	do{
			
		iter++;
			
		steps++;

		if (iter%100==0&&iter>0) {
			a_rate=(acceptpar+0.)/(trypar+0.);
			changepar=changepar*(0.95+a_rate);
			acceptpar=0;
			trypar=0;
		}

		for (int i=0; i<nvar; i++){

			move=(gsl_rng_uniform(rgen)*changepar)-(changepar/2);

			direction=gsl_rng_get(rgen);

			direction=(direction>0) ? 1:-1;
				
			move=direction*gsl_rng_uniform(rgen)*dxmax;
	
           		gsl_vector_set( xcurr, i,gsl_vector_get(xlastbest,i)+move);


           		lower=gsl_matrix_get(range,i,0);

           		upper=gsl_matrix_get(range,i,1);

           		if(gsl_vector_get(xcurr,i)<lower) gsl_vector_set(xcurr, i,lower+fabs(gsl_vector_get(xcurr,i)-lower));


           		if(gsl_vector_get(xcurr,i)>upper) gsl_vector_set(xcurr, i,upper-fabs(gsl_vector_get(xcurr,i)-upper));

			
		}

		trypar++;

		f2=obj_fn(xcurr,params); // compute likelihood for current parameters

		if(f2<f1){

			f1=f2;
			gsl_vector_memcpy(xlastbest,xcurr);
			acceptpar++;
		}

			
		// Convergence test

	         if (iter > min_iter && f1 != f2){
	                if (f1-f2 > 0.0 && f1-f2 < prec){
	                    ct++;
	                    if (ct==min_iter){ status = 0;}

			    f1=f2;
			    gsl_vector_memcpy(xlastbest,xcurr);
	                }
	                else{
	                    ct=0;
	                }

	         } else if(iter > min_iter && f1 == f2){
	            	ct++;

	                if (ct==min_iter){
                           status = 0;	
                        }

	                f1=f2;
			gsl_vector_memcpy(xlastbest,xcurr);
	         }
	}while (status &&  iter < max_iter);


	//Sammple around the maximum in order to test robustness of estimate

		
	for(int tt=0;tt<num_iter_pre;tt++){
					
		for(int i=0;i<nvar;i++){
								
			lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
			upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
								
			if(lower<upper){
				lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
				upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
				gsl_vector_set(xcurr_max,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
								
			}else{
				lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
				upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
				gsl_vector_set(xcurr_max,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
								
			}
					
		}

		f2=obj_fn(xcurr_max,params);


		if(f2<f1){
			f1=f2;
			gsl_vector_memcpy(xlastbest,xcurr_max);
		}
	}


// Store (update) best values ...
	

	for(int i=0;i<nvar;i++){
		other->data[i]=gsl_vector_get(xlastbest,i);
	}
	
// Cleanup...
	gsl_vector_free(xcurr);
	gsl_vector_free(xcurr_max);
	gsl_vector_free(xlastbest);

	return(f1);

}


// The  numerical minimization routine to get parameters by a hierarchical MC method...(only used for likelihood optimization under a K beta delay-deterministic model)...

double find_local_optimum_MChier (gsl_rng *rgen, gsl_vector * other, gsl_matrix * range,
		  void * params,
		  double (*obj_fn)( const gsl_vector * x, void * p),
		  double prec,
		  int& steps,
		  int verbose){

	int iter,max_iter =10000;
	int num_iter_pre=100000;
	int min_iter=1000;
	int status,ct=0;
	double f1=0, f2=0;

	//int direction;



	double move;

	int nvar = (int) other->size;
	gsl_vector *  xcurr=gsl_vector_alloc(nvar);//allocate here
	gsl_vector *  xcurr_max=gsl_vector_alloc(nvar);//allocate here
	gsl_vector *  xlastbest=gsl_vector_alloc(nvar);//allocate here

	//gsl_vector *  xbest = gsl_vector_alloc(nvar);

	//fpar * myfpar = static_cast<fpar*> (params);

	double upper=0.0;
	double lower=0.0;

	//double ftest=0.0;

	double dxmax=0.05;

	double changepar=dxmax;

	double a_rate=0.0;

	int acceptpar=0;

	int trypar=0;

	double rrr=0.0;

	status=1;

	iter=0;

	steps=0;

	ct=0;

	f2=1.0e20;

	// First fit a 1 beta delay-deterministic model //


	for(int tt=0;tt<num_iter_pre;tt++){

		for(int i=0;i<nvar;i++){


			lower=gsl_matrix_get(range,i,0);

			upper=gsl_matrix_get(range,i,1);

			if(i==0 || i==1){
				rrr=gsl_rng_uniform(rgen);
				gsl_vector_set(xcurr,i,lower+(rrr*(upper-lower)));
			}else{
				gsl_vector_set(xcurr,i,gsl_vector_get(xcurr,1));
			}
		}

		if(tt==0) {
			gsl_vector_memcpy(xlastbest,xcurr);
			f1=obj_fn(xcurr,params); // compute likelihood for current parameters
		}else{
			f2=obj_fn(xcurr,params);
		}

		if(f2<f1){
			f1=f2;
			gsl_vector_memcpy(xlastbest,xcurr);
		}

	}


		do{
			iter++;
			steps++;


			if (iter%100==0&&iter>0) {
				a_rate=(acceptpar+0.)/(trypar+0.);
				changepar=changepar*(0.95+a_rate);
				acceptpar=0;
				trypar=0;
			}

			for (int i=0; i<nvar; i++){

				move=(gsl_rng_uniform(rgen)*changepar)-(changepar/2);
				
		           	gsl_vector_set( xcurr, i,gsl_vector_get(xlastbest,i)+move);


           			lower=gsl_matrix_get(range,i,0);
           			upper=gsl_matrix_get(range,i,1);

           			if(gsl_vector_get(xcurr,i)<lower) gsl_vector_set(xcurr, i,lower+fabs(gsl_vector_get(xcurr,i)-lower));


           			if(gsl_vector_get(xcurr,i)>upper) gsl_vector_set(xcurr, i,upper-fabs(gsl_vector_get(xcurr,i)-upper));
			}

			trypar++;

			f2=obj_fn(xcurr,params); // compute likelihood for current parameters

			if(f2<f1){

				f1=f2;
				gsl_vector_memcpy(xlastbest,xcurr);


				acceptpar++;
			}

			// Convergence test

	         	if (iter > min_iter && f1 != f2){
	               	 	if (f1-f2 > 0.0 && f1-f2 < prec){
	                    		ct++;

	                    		if (ct==min_iter) status = 0;

			    		f1=f2;
			    		gsl_vector_memcpy(xlastbest,xcurr);

	                	}else{
	                    	ct=0;
	                	}

	         	} else if(iter > min_iter && f1 == f2){
	            		ct++;

	                	if (ct==min_iter){
                           	status = 0;	
                        	}

	                	f1=f2;
				gsl_vector_memcpy(xlastbest,xcurr);
	        	 }
		}while (status &&  iter < max_iter);




		//Sammple around the maximum in order to test robustness of estimate

		for(int tt=0;tt<num_iter_pre;tt++){
			for(int i=0;i<nvar;i++){
				lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
				upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
				if(lower<upper){
					lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
					upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
					gsl_vector_set(xcurr_max,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
				}else{
					lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
					upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
					gsl_vector_set(xcurr_max,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
				}
			}

			f2=obj_fn(xcurr_max,params);


			if(f2<f1){
				f1=f2;
				gsl_vector_memcpy(xlastbest,xcurr_max);
			}
		}


	// Then fit a K Beta //


		gsl_vector_memcpy(xcurr,xlastbest);

		for(int tt=0;tt<num_iter_pre;tt++){

			for(int i=0;i<nvar;i++){


				lower=gsl_vector_get(xcurr,i)-0.95*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));

				upper=gsl_vector_get(xcurr,i)+0.95*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));

				if(lower<upper){

					lower=gsl_vector_get(xcurr,i)-0.95*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
							
					if(tt<=90000){
						upper=gsl_vector_get(xcurr,i)+0.95*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
					}else{
						upper=gsl_matrix_get(range,i,1);
					}

					gsl_vector_set(xcurr,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));

				}else{
					if(tt<=90000){
						lower=gsl_vector_get(xcurr,i)-0.95*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
					}else{
						lower=gsl_matrix_get(range,i,0);
					}

					upper=gsl_vector_get(xcurr,i)+0.95*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));

					gsl_vector_set(xcurr,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));

				}

			}


			f2=obj_fn(xcurr,params);


			if(f2<f1){
				f1=f2;
				gsl_vector_memcpy(xlastbest,xcurr);
			}

		}


				do{
					iter++;
					steps++;
		
		
					if (iter%100==0&&iter>0) {
						double a_rate=(acceptpar+0.)/(trypar+0.);
						changepar=changepar*(0.95+a_rate);
						acceptpar=0;
						trypar=0;
					}
		
					for (int i=0; i<nvar; i++){
		
						move=(gsl_rng_uniform(rgen)*changepar)-(changepar/2);
			
		           			gsl_vector_set( xcurr, i,gsl_vector_get(xlastbest,i)+move);
		
			           		lower=gsl_matrix_get(range,i,0);
		
			           		upper=gsl_matrix_get(range,i,1);
		
			           		if(gsl_vector_get(xcurr,i)<lower) gsl_vector_set(xcurr, i,lower+fabs(gsl_vector_get(xcurr,i)-lower));
		
		
			           		if(gsl_vector_get(xcurr,i)>upper) gsl_vector_set(xcurr, i,upper-fabs(gsl_vector_get(xcurr,i)-upper));
		
					}
		
					trypar++;
		
		
		
					f2=obj_fn(xcurr,params); // compute likelihood for current parameters
		
					if(f2<f1){
		
						f1=f2;
						gsl_vector_memcpy(xlastbest,xcurr);
		
						acceptpar++;
					}
					// Convergence test

	         			if (iter > min_iter && f1 != f2){
	               	 			if (f1-f2 > 0.0 && f1-f2 < prec){
	                    				ct++;

	                    				if (ct==min_iter){ status = 0;}

			    				f1=f2;
			    				gsl_vector_memcpy(xlastbest,xcurr);

	                			}else{
	                    				ct=0;
	                			}

	         			}else if(iter > min_iter && f1 == f2){
	            					ct++;

	                				if (ct==min_iter){
                           				status = 0;	
                        				}

	                				f1=f2;
							gsl_vector_memcpy(xlastbest,xcurr);
	        	 		}

				}while (status &&  iter < max_iter);




				//Sample around the maximum
				for(int tt=0;tt<num_iter_pre;tt++){
					for(int i=0;i<nvar;i++){
						lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
						upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));

						if(lower<upper){
							lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
							upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,0));
							gsl_vector_set(xcurr_max,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
						}else{
							lower=gsl_vector_get(xcurr,i)-0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
							upper=gsl_vector_get(xcurr,i)+0.25*fabs(gsl_vector_get(xcurr,i)-gsl_matrix_get(range,i,1));
							gsl_vector_set(xcurr_max,i,lower+(gsl_rng_uniform(rgen)*(upper-lower)));
										
						}

						f2=obj_fn(xcurr_max,params);


						if(f2<f1){
							f1=f2;
							gsl_vector_memcpy(xlastbest,xcurr_max);
						}
					}
				}



	// Store (update) best values ...


	for(int i=0;i<nvar;i++){
		other->data[i]=gsl_vector_get(xlastbest,i);
	}

// cleanup...

	gsl_vector_free(xcurr);
	gsl_vector_free(xcurr_max);
	gsl_vector_free(xlastbest);

	return(f1);

}



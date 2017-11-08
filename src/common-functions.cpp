/*

 common-functions.cpp
 
 Created on: Nov 2017
  
 Author: Nuno Rocha Nene
 E-mail:  nunonene@gmail.com
 
*/

#include "common-functions.h"
#include "emission.h"


using namespace std;


// Get general dimensions of a data set for...

int get_dims(const char * data_fn,int& maxNsim, int& Khaplo,int& T){

  ifstream data_ifs;

  string line;

  stringstream line_ss;

  ostringstream convert;   // stream used for the conversion

  string h="_Haplo";

  string smpl="_smpl";

  string end;

  end=".txt";

  string file;
 
  convert << (Khaplo-1); // append string; 
	
  string hnumb=convert.str();

  convert.str("");
 
  convert.clear();

  file=data_fn+h+hnumb+smpl+end; // append string



  //Check if file oppened

  data_ifs.open( file.c_str(), ios::in);

  if (data_ifs.fail()){
  		    
	printf("ERROR: file %s cannot be opened.\n", file.c_str());

        exit(1);  	
  }


  double r;

  int nT=0; // Number of generations; Needs to be equal to T;

  int nSim=0; // Number of lines in each files; Needs to be equal to paramter Nsim provided at the command line;

  int count1=0;
  int count2=0;
  int count=0;

  while( data_ifs.good()){

    line.clear();
    getline( data_ifs, line);

    if (line.empty()) break;

    line_ss.clear();

    line_ss.str(line);

    nSim++;

    //Check if the number of time-points per line, i.e. per simulation, is consistent with value provided at the command line;

    nT=0;

    count1=0;

      while(line_ss >> r){
    	  nT++;

    	  if(r==1)count1++;
    	  
  	  if(count1==50)count2=nT;
      }

      if(T!=nT){
      	 // printf("ERROR: number of  generations in %s, line %i, is not consistent with total duration T .\n", data_fn,nSim);
	    printf("ERROR: number of  generations in %s, line %i, is not consistent with total duration T .\n", file.c_str(),nSim);
            exit(1);
        }

      if(count2>count)count=count2;
  }

  data_ifs.close();

  if(nSim<maxNsim){
         	 // printf("ERROR: number of  simulations in %s is smaller than maxNsim=%i .\n", data_fn,maxNsim);
 		  printf("ERROR: number of  simulations in %s is smaller than maxNsim=%i .\n", file.c_str(),maxNsim);
               exit(1);
   }
   data_ifs.close();

   if(count==0)count=T;

   return(count);
}


// Read in data: Nsim rows and T columns; each file corresponds to a different haplotype

void get_data( const char * data_fn_smpl, Emission * myEmit, vector<int>& selectlines){

  ifstream data_ifs;
  string line;
  stringstream line_ss;

  ostringstream convert;

  string h="_Haplo";

  string smpl="_smpl";
  string end;

  end=".txt";

  
  string hnumb;

  string file;

       
  for (int i=0;i<myEmit->maxNsim;i++){ // at the moment selects maxNsim number of trajectories has they appear in file

    	   selectlines.push_back(i);
  }

  for (int k=0;k<myEmit->KH;k++){
	 
	convert.clear();

	convert<<k;

	hnumb=convert.str();

	convert.str("");
	  	  
	convert.clear();
	  
	file=data_fn_smpl+h+hnumb+smpl+end; // append string

	data_ifs.open(file.c_str(), ios::in); // open file

	// Check if file is opened

	if (data_ifs.fail()){
		  
	printf("ERROR: file %s cannot be opened.\n", file.c_str());
		 
	exit(1);
	  
	}

	  
	int ct=0;
	  
	int ns=0;
	  
	int sel=0;
	  
	int st=0;
	  
	double r;


 
    	// Collect all data...

	 while( data_ifs.good()){		
		
		line.clear();
		  
		getline( data_ifs, line);
		  
		if (line.empty()) break;
		  
		line_ss.clear();
		  
		line_ss.str(line);
 
		if(sel<myEmit->maxNsim){
			if(ns==selectlines[sel]){				  
				ct=0;
				st=0;

				while(line_ss >> r){

					  if (st==0){ // read only data at the sampling genrations

						  myEmit->reads[k][sel][st] = floor((double)r*(double)(myEmit->C));

						  st++;

					  } else if (st< myEmit->nTs){

						  if (ct == myEmit->sTimes[st]) {						 
							  
							myEmit->reads[k][sel][st] = floor((double)r*(double)(myEmit->C));

							st++;
						  }

					  }
					  ct++;
				  }
				  sel++;
			
			  }

		  }
		  ns++; // count one simulation read
	  }
    
		data_ifs.close();

  }

}

// Get weighted mutation graph for the linear fitness chain or the fitness hypercube

void get_mutgraph(gsl_matrix * MutProp,gsl_vector * FSPL,int Type,int KH){



    if(Type==0){ // if linear chain

    		double val=0.0;

    		for (int k1=0; k1<KH; k1++){
    			for (int k2=0; k2<KH; k2++){
    				val=0.0;
    				if(k1==k2) {
    					val=-1;
    				}else{
    					if(k1==1 && k2==0){
    						val=1;
    					}
    					if(k1==KH-2 && k2==KH-1){
    						val=1;
    					}
    					if(KH>2){
    						if(k1==1 && k2==k1+1){
    							val=0.5;
    						}

    						if(k1==(KH-1) && k2==(k1-1)){
    							val=0.5;
    						}

    						if(k1!=1 && k1!=(KH-2) && (k2==k1+1 || k2==k1-1)){
    							val=0.5;
    						}

    						if(k1==(KH-2) && k2==k1-1){
    							val=0.5;
    						}
    					}
    				}
    				gsl_matrix_set(MutProp,k1,k2,val);
    			}
    		}
    }else{

    // If hypercube

       get_HyperGraphMutMat(MutProp,KH); // Get mutation matrix weights

    

       get_HyperGraphSPL(FSPL,KH,Type); // Get shortest path distance on graph

      

    }





}



// Get hypergraph mutation matrix

void get_HyperGraphMutMat(gsl_matrix * MutProp,int n){

	ifstream data_ifs;

	string line;

	stringstream line_ss;
	 
	ostringstream convert;
	  
	string data_fn="../HyperGraphMutFit/HyperGraphMutMat";
 
	string h="_nh";
	string l=".txt";

	string hnumb;

	string file;

	convert << (n);

	hnumb=convert.str();

	file=data_fn+h+hnumb+l; // append string
		  
	convert.str("");
		  
	convert.clear();

        data_ifs.open(file.c_str(), ios::in); // open file
	  
	// Check if file is opened ...
	if (data_ifs.fail()){
	  
	  printf("ERROR: file %s cannot be opened.\n", file.c_str());
	  exit(1);
		  
	}

		 
	int r;
	int c;
	double val;

	// Now collect all data...

		  
	while( data_ifs.good()){
			  
		line.clear();
			  
		getline( data_ifs, line);
			  
		if (line.empty()) break;
			  
		line_ss.clear();
			  
		line_ss.str(line);
	 
		while(line_ss >> r >> c >> val){

		gsl_matrix_set(MutProp,r-1,c-1,val);
			  
		}	  
	}

	data_ifs.close();

}

// Collect distance to peak (Hypercube only)

void get_HyperGraphSPL(gsl_vector * FSPL,int n,int Type){

	ifstream data_ifs;

	string line;

	stringstream line_ss;

	ostringstream convert;
	  	  
	string data_fn="../HyperGraphMutFit/HyperGraphMutMatSPL";
	
	  
	string h1="_nh";

	string l=".txt";

	string hnumb1;

	string file;

	convert << (n);

	hnumb1=convert.str();

	convert.str("");
		  
	convert.clear();

	file=data_fn+h1+hnumb1+l; // append string

	data_ifs.open(file.c_str(), ios::in); // open file

	if (data_ifs.fail()){
		printf("ERROR: file %s cannot be opened.\n", file.c_str());
		exit(1);
	}

	double val;
	int count=0;

	// Now collect all data...

       while( data_ifs.good()){
		line.clear();
		getline( data_ifs, line);
			  
		if (line.empty()) break;
			  
		line_ss.clear();
			  
		line_ss.str(line);
  
		while(line_ss >> val){

		gsl_vector_set(FSPL,count,val);
		}

		count++;
	}
	    
	data_ifs.close();

}


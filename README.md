# Inference-of-fitness-effects-from-time-resolved-genome-sequence-data


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
 
 4. In order to compile see Makefile (requires gsl) in src folder. Once compiled copy ./FKHoptm to the Examples folder in order to run the code with the files provided.

 Example (see data in Examples folder):

 (Stochastic model)

 ./FKHoptm --data SimKhaploDynamics_Khaplo6_N1000000000_Nsim100_C100_Init1_sigma0.1_mu0.00001_T2000_TypeLin --out  SimKhaploDynamics_Khaplo6_N1000000000_Nsim100_C100_Init1_sigma0.1_mu0.00001_T2000_TypeLin --N 1000000000 --C 100 --mu 0.00001 --Nsim 100 --maxNsim 100 --T 2000 --nTs 200 --Khaplo 6 --si 0.01 --sf 1 --sgrid 100 --DSModel 2 --mode 1 --optmode 1 --Type 0;

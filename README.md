# pricing
Smoothly state dependent pricing

Program files for "Distributional Dynamics under Smoothly State-Dependent Pricing", 
   Bank of Spain working paper 0831, December 2008, (C) James Costain and Anton Nakov


OPERATING PROCEDURE 


1. Run "gedyn.m" to compute general equilibrium dynamics; Impulse-responses are computed automatically

2. Run "vd.m" to compute the variance decomposition and estimate Phillps curve coefficients

The default setting is a Taylor rule. To switch to a money growth rule, set "phiPI" to 0 in "gedyn".

The default model is SSDP. To change the model, change "adjtype" in "gedyn".


COMPUTATIONAL DETAILS

General equilibrium steady state calculation requires a few seconds on an ordinary Pentium 4 with 3Ghz CPU and 1GB of memory. 

Calculation of the dynamics requires approximately 2-3 minutes. This is hard to speed up since the bulk of the time (the slowest step) is spent on the QZ decomposition, which uses the in-built MATLAB function qz.


LIST OF ALL PROGRAM FILES


   adjustment        - Computes adjustment probability as a function of the gain from adjustment
   calcstats         - Computes steady-state statistics 
   compute_IRFs      - Computes dynamic paths as a function of initial state and shocks (Taylor rule)
   compute_IRFsM     - Computes dynamic paths as function of initial state and shocks (money rule)
   disclyap          - Solves discrete Lyapunov equation
   distsim           - Computes impulse-responses based on Klein's state-space solution
   dyneqklein        - Defines equations of stochastic model (Taylor rule)
   dyneqkleinM       - Defines equations of stochastic model (money rule)
   dynsolveklein     - Solves dynamic general equilibrium (Taylor rule)
   dynsolvekleinM    - Solves dynamic general equilibrium (money rule)
   estparam          - Sets idiosyncratic shock and adjustment function parameters 
   ExpectMenuCost    - Expected menu cost in SMC model
   findPEVfun        - Finds partial equilibrium steady-state value function for a given wage 
   gedyn             - Main program for computing general equilibrium dynamics 
   hazard            - Computes the adjustment hazard function after an initial price change
   histpchanges      - Computes the histogram of price changes
   histsimaggidio    - Simulates price histories with idiosyncratic productivity shocks 
   hpfilter          - Extracts the Hodrick-Prescott trend of a time series
   infdecomp         - Decomposes inflation into intensive, extensive, selection effects
   infparam          - Loads alternative gross money growth rates
   irf               - Computes and plots impulse-response functions
   jacob_reiter      - Computes Jacobian by forward differences
   kleinsolve        - Implements Klein's QZ decomposition method for solving linear RE models
   ks                - Computes Kolmogorov-Smirnov statistic for equality of two cdf's
   lamcontin         - Approximates lambda in the menu cost model with a continuous function 
   M_pStar           - Computes value function maximum and argmax with quadratic interpolation on V
   makegrids         - Builds productivity and price grids in logs       
   optimprint        - Prints current parameter vector at each iteration of the estimation procedure
   param             - Sets program execution and macro model parameters
   Pdist_iter        - Computes the stationary distributions of firms before and after shocks
   Pidentity         - Function called by fzero to find the general equilibrium steady-state
   plot_IRFs         - Plots impulse-response functions
   plotfigs          - Plots figures of steady-state objects
   Pmatrix           - Calculates matrix P which rounds stochastically to two grid points around the optimal price  
   printstats        - Prints out steady-state statistics
   progreport        - Reports convergence statistics
   Rmatrix           - Calculates the transition matrix R 
   setmat            - Initializes matrices based on guess for wage=wflex 
   sizepchtime       - Computes the evolution of the mean absolute price change after an initial adjustment
   tauchen           - Converts a VAR(1) into a Markov-Chain using Tauchen's method
   taylor            - Taylor rule
   V_iter            - Solves value function by iteration on a grid with interpolation
   vd                - Computes variance decomposition and Phillips curve regression


We thank Elmar Mertens for the code that implements Tauchen's method for approximating AR1 by a finite-state Markov process available on his website http://www.elmarmertens.ch/


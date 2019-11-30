% Sets program execution and macro model parameters
global mu

% PROGRAM EXECUTION PARAMETERS               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  idioshocks   = 1;                          % heterogeneity: 1-idio shocks; 0-fixed heterogeneity; -1 rep agent
  showconverge = 1*(version<3);              % show convergence progress. Set to 0 when estimating (version 3)
  
% DISCOUNT AND PREFERENCE PARAMETERS         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  beta    = (1.04^(-1/12));                  % utility discount factor  (raise to ^(7/30) for weekly frequency)
  gamma   = 2;                               % CRRA coefficient
  epsilon = 7;                               % elasticity of substitution among differentiated goods
  chi     = 6;                               % labor supply parameter
  nu      = 1;                               % money demand parameter
 
% STEADY STATE GROWTH RATE OF MONEY  
  mu      = infparam('1');                   % long-run gross inflation target   
  
% REAL WAGE GUESS  
  wflex   = (epsilon-1)/epsilon;             % guess for real wage based on flex price representative agent model
    
% CONVERGENCE TOLERANCE LEVELS               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ptol = (10^(8-accuracy))*eps;              % fzero convergence tolerance for aggregate price level identity
  Vtol = (10^(8-accuracy))*eps;              % convergence tolerance for Value function iteration
  PdistDIFFtol = (10^(8-accuracy))*eps;      % convergence tolerance for Pdist function iteration

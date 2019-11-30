% Main program for computing general equilibrium dynamics 
% Programs for "Dynamics of the Price Distribution in a General Model of
% State-Dependent Pricing", Working paper 0831, Bank of Spain
% Copyright (C) James Costain and Anton Nakov (2008)

global fzeroiter

% PROGRAM EXECUTION PARAMETERS
adjtype    = 0;       % 0-SSDP 2-Calvo 3-Fixed menu cost 4-Woodford 5-SMC
version    = 1;       % 1 model with estimated productivity 
finegrid   = 0;       % choose fine (1) or coarse (0) grid
gridSpread = 0.15;    % extra spread for price grid as share of (PMAX-PMIN)
accuracy   = 1;       % SS accuracy from 0 (lowest) to 4 (highest, slowest)
fzeroiter  = 0;       % reset iteration counter

% AGGREGATE SHOCK PARAMETERS  
phiM    = 0;          % persistence of monetary shocks
phiA    = 0.95;       % persistence of aggregate technology shocks 

% TAYLOR RULE PARAMETERS      
phiR    = 0.9;        % interest rate smoothing coefficient 
phiPI   = 2;          % inflation response coefficient
phiC    = 0.5;        % output response coefficient 

param;                % load macro parameters and guess for wages "wflex"
estparam;             % load the adj. function and productivity parameters

paramvec = [lbar, alpha, ksi, rho, stdMC, lam0]';  % collect in vector

options = optimset('Display','off','TolFun',ptol); % set fzero options 

% find equilibrium real wage with "wflex" as the inital guess
[wbar,fval,exitflag] = fzero('Pidentity',wflex,options,paramvec,version,...
                                adjtype,gridSpread,finegrid,ptol,accuracy);

if abs(fval)>ptol, error('fzero convergence failure'); end

load GE_ss;          % load steady-state into workspace
calcstats;           % calculate steady-state statistics

% compute dynamics using Klein's method
if phiPI>0, dynsolveklein    % Taylor rule 
else        dynsolvekleinM   % money growth rule
end

irf;                 % plot impulse-response functions
delete GE_ss.mat;    % delete file to save disk space
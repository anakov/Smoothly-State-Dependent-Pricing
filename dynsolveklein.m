% Solves dynamic general equilibrium (Taylor rule)
% dynsolveklein.m calls dyneqklein.m, jacob_reiter.m, and kleinsolv.m.

global Params;
format long
jacstep = 100*ptol;                   
eqcutoff = 100*jacstep;

phiMAT   = [phiM 0; 0 phiA];       % Matrix of shock persistence parameters

% X =[Pdist, DP, Rlag, V, C, PI]
nz  = 2;                 % number of exogenous shock processes (money, TFP)
nsj = 3;                 % number of scalar jump variables (C, PI, R)
nss = 2;                 % number of scalar state variables (DP, Rlag)
nPhi = gridsize;         % number of points of the distribution (states)
nV   = gridsize;         % number of points of the value function (jumps)
Ntot = nPhi+nss+nV+nsj;  % total number of variables            

% STORE VARIABLES IN PARAMS
Params.nz = nz;
Params.nsj = nsj;
Params.nss = nss;
Params.nPhi = nPhi;
Params.nV = nV; 

% Y WILL CONSIST OF FOUR PARTS
Params.ix     = 1:Ntot;                   % length Ntot
Params.ixnext = Ntot+1:2*Ntot;            % length Ntot
Params.iz     = 2*Ntot+1:2*Ntot+nz;       % length nz
Params.iznext = 2*Ntot+nz+1:2*Ntot+2*nz;  % length nz

Params.beta = beta; 
Params.gamma = gamma; 
Params.chi = chi; 
Params.epsilon= epsilon; 
Params.mu = mu; 
Params.nu = nu; 

Params.Rss = Rss; 
Params.Cbar = Cbar; 
Params.phiPI = phiPI; 
Params.phiC = phiC; 
Params.phiR = phiR; 

Params.adjtype = adjtype; 
Params.lbar = lbar; 
Params.lam0 = lam0; 
Params.ksi = ksi;              
Params.alpha = alpha;          

Params.PMAT = PMAT; 
Params.sMAT = sMAT; 
Params.TRANSMAT = TRANSMAT; 
Params.Pgrid = Pgrid; 
Params.pstep = pstep; 

% X = [Pdist, DP, Rlag, V, C, PI]
Xss = [Pdist(:);WeightPriceDispers;Rss;V(:);Cbar;mu;Rss];   
stst2=[Xss;Xss;zeros(nz,1);zeros(nz,1)];

resid = feval(@dyneqklein,stst2);
check = max(abs(resid));
   if(check>eps^.5 || 10*check>jacstep)
     disp('WARNING: Large residual at stst:');
     disp(resid);
     keyboard
   end;
tic
disp(sprintf('\n'))  
disp('START JACOBIAN CALCULATION')

   jac = jacob_reiter(@dyneqklein,stst2,jacstep);
   ajac = jac(:,Params.ixnext);
   bjac = jac(:,Params.ix);
   cjac = jac(:,Params.iznext);
   djac = jac(:,Params.iz);

disp(sprintf('\n'))  
disp('START KLEIN SOLUTION')

[JUMPS,STATEDYNAMICS,stableeigs,unstableeigs] = ...
    kleinsolve(ajac,bjac,cjac,djac,phiMAT,nPhi+nss,eqcutoff);

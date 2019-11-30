% Defines equations of stochastic model (Taylor rule)
function Resid = dyneqklein(Y)
global Params;
% Input: X: All variables, current and lagged
% Outputs: Equation residuals
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD parameters
 
beta = Params.beta; 
gamma = Params.gamma; 
chi = Params.chi; 
epsilon= Params.epsilon; 
mu = Params.mu; 
% nu = Params.nu; 

Rss = Params.Rss; 
Cbar = Params.Cbar; 
phiPI = Params.phiPI; 
phiC = Params.phiC; 
phiR = Params.phiR; 

adjtype = Params.adjtype; 
alpha = Params.alpha; 
lbar = Params.lbar; 
lam0 = Params.lam0; 
ksi = Params.ksi; 

PMAT = Params.PMAT; 
sMAT = Params.sMAT; 
TRANSMAT = Params.TRANSMAT; 
Pgrid = Params.Pgrid; 
pstep = Params.pstep; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KLEIN setup: X = [Pdist, DP, V, C, PI, Rlag]
% z is the exogenous stochastic processes driving money and productivity, with z_t+1 = phiz z_t + iidshock
% PhiHatNow is dist at beginning of t before firms adjust AND before money shock z_t is realized.
% PhiHatNow_t = Phi_t-1 = end of period t-1 production distribution
% State variables are: today's shocks z, price dispersion and PhiHatNow  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD VARIABLES 
  nPhi = Params.nPhi;
  nV = Params.nV;
% nsj = Params.nsj;
  nss = Params.nss;

% EXOGENOUS VARIABLES  
  znow = Y(Params.iz);       % exogenous shock process (possibly correlated)
  znext = Y(Params.iznext);  % ShockNow and ShockNext are scalars: no further extraction needed
  zRnow = znow(1);           % Taylor rule shock today
% zRnext = znext(1);         % Taylor rule shock tomorrow
  zAnow = znow(2);           % Common TFP shock

% ENDOGENOUS VARIABLES:
  xnow = Y(Params.ix);         % endogenous variables now
  xnext = Y(Params.ixnext);    % endogenous variables next
    
% ENDOGENOUS STATE VARIABLES
  PhiHatNow = xnow(1:nPhi);
  PhiHatNext = xnext(1:nPhi);
% DPnow = xnow(nPhi+1);
  DPnext = xnext(nPhi+1);
  ilagnow = xnow(nPhi+2);  % lag of nominal interest rate
  ilagnext = xnext(nPhi+2);

  % ENDOGENOUS JUMP VARIABLES:
  Vnow = xnow(nPhi+nss+1:nPhi+nss+nV);
  Vnext = xnext(nPhi+nss+1:nPhi+nss+nV);
  Cnow = xnow(nPhi+nss+nV+1);
  Cnext = xnext(nPhi+nss+nV+1);
  PInow = xnow(nPhi+nss+nV+2);
  PInext = xnext(nPhi+nss+nV+2);
  inow = xnow(nPhi+nss+nV+3);
  inext = xnext(nPhi+nss+nV+3);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESHAPE VECTORS TO MATRICES size (nump,nums)
  [nump,nums]=size(PMAT);
  siz=[nump,nums];
  PhiHatNow = reshape(PhiHatNow,siz);
  PhiHatNext = reshape(PhiHatNext,siz);
  Vnow = reshape(Vnow,siz);
  Vnext = reshape(Vnext,siz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE SOME VARIABLES USING EQUATIONS OUTSIDE OF KLEIN REPRESENTATION
  
% Real wage: from Labor supply equation (detrended by P)
  wnow  = chi*Cnow^gamma; 
  wnext = chi*Cnext^gamma;
  
% Payoff today
  Unow = Cnow*(PMAT).^(-epsilon).*(PMAT - wnow*exp(-zAnow).*sMAT);   
  
% Calculate the transition matrices R (THE OFFSET IS DETERMINED BY PI NOW)
  Rnow  = Rmatrix(nump, PInow, pstep);
  Rnext = Rmatrix(nump, PInext, pstep);
 
% Get optimal value and optimal price now and next
  [Mnow, pStarNow]   = M_pStar(Vnow, Pgrid, PInow);
  [Mnext, pStarNext] = M_pStar(Vnext, Pgrid, PInext);

% Calculate adjustment values Dnow and Dnext
  Dnow  = ones(nump,1)*Mnow - Vnow;   
  Dnext = ones(nump,1)*Mnext - Vnext; 
  
% Calculate probabilities of adjustment Lambdanow and Lambdanext, and ExpectGAINSnext    
  LambdaNow = adjustment(adjtype, Vnow, Dnow, wnow, ksi, alpha, lbar, lam0);
  [LambdaNext ExpectGAINSnext] = adjustment(adjtype, Vnext, Dnext, wnext, ksi, alpha, lbar, lam0);

% Calculate PhiTildeNow from PhiHatNow = Phi_t-1
  PhiTildeNow = Rnow*PhiHatNow*TRANSMAT';
  
% Calculate matrix P indicating the split of adjusting firms' 
% mass between adjacent grid points around the optimal price  
  PmatNow = Pmatrix(pStarNow, Pgrid, pstep);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE RESIDUALS FOR JACOBIAN

% PhiHatNext = dist of prices at time of production NOW (after shocks and adjustment)
  PhiResid = -PhiHatNext + (1-LambdaNow).*PhiTildeNow + PmatNow.*(ones(nump,nump)*(LambdaNow.*PhiTildeNow)); 

% Value function residual
  VResid = -Vnow + Unow + beta*(Cnext/Cnow)^(-gamma)*Rnext'*ExpectGAINSnext*TRANSMAT;
 
% Consumption Euler residual
  eulerResid = 1 - inow*beta*(Cnext/Cnow)^(-gamma)/PInext; 

% Aggregate price level residual
  priceResid = 1 - (sum(sum(PhiHatNext.*(PMAT.^(1-epsilon)))))^(1/(1-epsilon));
  
% Price dispersion dynamics 
  deltaResid = DPnext - sum(sum(PhiHatNext.*sMAT.*PMAT.^(-epsilon))); 
   
% Lagged interest rate identity
  ilagresid = ilagnext - inow;
  
% Nominal interest rate: from Taylor rule
  iresid = inow - taylor(Rss,phiR,ilagnow,phiPI,PInow,mu,phiC,Cnow,Cbar,zRnow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE RESID VECTOR:
  Resid = [PhiResid(:); VResid(:); eulerResid; priceResid; deltaResid; ilagresid; iresid];
  % order of residuals is unimportant, but for simplicity follow same ordering
  % for residuals as for variables and equations
   

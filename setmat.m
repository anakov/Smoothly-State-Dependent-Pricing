% Initializes matrices based on guess for wage=wflex 

% Matrices representing idiosyncratic state: converted to LEVELS
sMAT=ones(nump,1)*exp(sgrid);  % exogenous (INVERSE) idiosyncratic productivity with transition matrix TRANSMAT
PMAT=exp(Pgrid)*ones(1,nums);  % endogenous sticky idiosyncratic price

% Implied value of C based on the guess for w
Cflex = (wflex/chi)^(1/gamma);                 

% Implied value of the real PAYOFF based on the guess for w and Cflex
PAYOFFflex = Cflex*PMAT.^(-epsilon).*(PMAT-wflex*sMAT);

% Initial guess for value function     
Vflex = PAYOFFflex/(1-beta);                      

% Initial guess for the distribution of firms    
Pdist=zeros(nump,nums);          

if idioshocks==0 % fixed heterogeneity
   % Pdist=ones(nump,nums)/gridsize;       % uniform distribution  
   sgridcount = normpdf(sgrid,0,stdMC);    % normally distributed firm
   sgridcount = sgridcount/sum(sgridcount);
   pgridcount = normpdf(Pgrid,0,stdMC); 
   pgridcount = pgridcount/sum(pgridcount);
   Pdist = pgridcount*sgridcount;
else
   Pdist(ceil(nump/2),ceil(nums/2))=1;   % initial probability distribution with unit atom in middle
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%BUILDING TRANSMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each ROW of TRANSMAT refers to some state at t+1, while each COL of TRANSMAT refers to some state at t.
sigma2_eps = stdMC^2*(1-rho^2);                 % variance of technology innovation  
if idioshocks==1
  [TRANSMAT, out1, out2, out3, out4, out5] = tauchen(rho,0,sigma2_eps,nums,numstddevs);

  TRANSMAT=TRANSMAT';                         % Productivity state transition matrix
  if any(abs(out1'-sgrid)>eps^.5), error('Grid mismatch'); end
  if any(abs(sum(TRANSMAT)-1)>eps^.5), error('A column of TRANSMAT does not sum to 1'); end

elseif idioshocks==0
  TRANSMAT=eye(nums);

elseif idioshocks==-1
  if rem(nums,2)==0, error('nums must be odd for rep agent case'), end    
  TRANSMAT=zeros(nums,nums);
  TRANSMAT((1+nums)/2,:)=ones(1,nums);
end
% COLUMNS of TRANSMAT should sum to one.

%BUILDING RECURSEMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RECURSEMAT = Rmatrix(nump, mu, pstep);  
% RECURSEMAT is defined so that its COLUMNS sum to one
% Each ROW of RECURSEMAT refers to some price at t+1, and 
% each COL of RECURSEMAT refers to some price at t.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


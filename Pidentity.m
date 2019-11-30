function zero = ... 
 Pidentity(wbar,paramvec,version,adjtype,gridSpread,finegrid,ptol,accuracy)
% Function called by fzero to find the general equilibrium steady-state
% searches over the wage rate that satisfies the aggregate price identity
% input : guess for wages
% output: residual of price identity equation (fzero minimizes this)

lbar  = paramvec(1);
alpha = paramvec(2);
ksi   = paramvec(3);
rho   = paramvec(4);
stdMC = paramvec(5);

if adjtype==8, lam0 = paramvec(6); else lam0 = []; end

findPEVfun;      % find steady-state value function for a given wage "wbar"

while (any(pstar == PMIN) || any(pstar == PMAX))  % if grid boundary is hit
  disp('Corner solution found. Hit any key to stretch the grid.')  
  pause
  gridSpread = gridSpread + 0.1;     % then stretch out the price grid
  findPEVfun;                        % find steady-state on extended grid
end

Pdist_iter;                          % get stationary distribution of firms

zero = 1-sum(sum(PMAT.^(1-epsilon).*Pdist))^(1/(1-epsilon)); % get residual

save GE_ss
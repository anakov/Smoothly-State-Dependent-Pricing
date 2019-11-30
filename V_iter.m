% Solves value function by iteration on a grid with interpolation
% version 17 april 2008
global fzeroiter

[nump,nums]=size(V);

VDIFF=inf;                                                 % reset difference in value function
Viter=0;                                                   % reset counter
fzeroiter=fzeroiter+1;

while VDIFF>Vtol && Viter<25000                            % iterate to convergence of V 
  Viter=Viter+1;                                           % increment counter
  if rem(Viter,100)==0  && showconverge                    % report convergence progress
     progreport(adjtype,finegrid,fzeroiter,VDIFF,Viter); 
  end

% Time t+1 values

  [M, pstar] = max(V);                                     % discrete maximization
  if  min(pstar)>1 && max(pstar)<nump                      % if the solution is interior
    [M, pstar] = M_pStar(V, Pgrid, mu);                    % quadratic value function interpolation at each state
  end
  
  D = ones(nump,1)*M - V;                                  % D is the value of adjustment

  if any(any(D<-eps^.5)),                                  % D should be non-negative
     % disp('Negative D after interpolation of V.')
      D(D<-eps^.5)=0;
  end

  [lambda, EXPGAINS] = adjustment(adjtype, V, D, wbar, ksi, alpha, lbar, lam0);

% Time t values
  iterV = PAYOFFMAT + beta*RECURSEMAT'*EXPGAINS*TRANSMAT;   % iterV is current payoff plus disc. continuation value
  VdifMin = min(min(iterV-V));                             % minimum difference 
  VdifMax = max(max(iterV-V));                             % maximum difference
  VDIFF = VdifMax-VdifMin;                                 % distance between max and min difference
 % VDIFF = max(max(abs(iterV-V)));                            % change in value function (sup norm)
  V = iterV;                                               % updating V
end
 V = V + beta/(1-beta)*(VdifMax+VdifMin)/2;                % vertical shift of value function once shape has converged

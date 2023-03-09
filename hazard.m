function hazd = hazard(horizon, RECURSEMAT, Pdist, TRANSMAT, P, lambda)
% Computes the adjustment hazard function after an initial price change
% hazd(k) = fracadjusters(k)/fracsurvivors(k-1) 
% conditional on all firms having adjusted at time 0

[nump, out] = size(Pdist);
fracsurvivors = NaN*ones(horizon+1,1);
hazd = NaN*ones(horizon+1,1);

% period 0, forcing price adjustment from steady-state
PDE = RECURSEMAT*Pdist*TRANSMAT';                              % distribution after productivity shock and inflation
AdjusterMAT = P.*(ones(nump,nump)*(lambda.*PDE));              % distribution of adjusters after adjustment
SurvivorMAT = AdjusterMAT;                                     % survivors at time 0 is everybody who adjusted at 0
fracsurvivors(1) = sum(sum(SurvivorMAT));                      % fraction of survivors at time 0

% periods from 1 onwards: keeping track of fracsurvivors only
for k=2:horizon+1
 PDE = RECURSEMAT*SurvivorMAT*TRANSMAT';                       % distribution after productivity shock and deflation
 SurvivorMAT = (1-lambda).*PDE;                                % distribution of fracsurvivors at k 
 fracadjusters  = sum(sum(lambda.*PDE));                       % fraction of adjusters at k 
 fracsurvivors(k) = sum(sum(SurvivorMAT));                     % fraction of survivors at k
 hazd(k) = fracadjusters/fracsurvivors(k-1);                   % hazard at k 
 end

hazd = hazd(2:end);                                            % eliminating time 0 value

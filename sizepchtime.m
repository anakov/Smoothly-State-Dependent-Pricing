function MeanAbsPchange = sizepchtime(horizon, AbsPdistans,RECURSEMAT,Pdist,P,TRANSMAT,lambda)
% Computes the evolution of the mean absolute price change after an initial adjustment

[nump, out] = size(Pdist);
MeanAbsPchange = NaN*ones(horizon+1,1);                             % preallocate memory for the vector

% period 0, forcing price adjustment from steady-state
PDE = RECURSEMAT*Pdist*TRANSMAT';                                   % distribution after shocks
AdjusterMAT = P.*(ones(nump,nump)*(lambda.*PDE));                   % distribution of adjusters after adjustment
SurvivorMAT = AdjusterMAT;                                          % survivors at 0 (everyone who adjusted at 0)

% periods from 1 onwards: keeping track of survivors only
for k=1:horizon+1
PDE = RECURSEMAT*SurvivorMAT*TRANSMAT';                             % distribution of firms after shocks
AdjusterMAT = lambda.*PDE;                                          % distribution of adjusting firms
fracadjusters = sum(sum(AdjusterMAT));                              % fraction of adjusters at time k
MeanAbsPchange(k)=sum(sum(AbsPdistans.*AdjusterMAT))/fracadjusters; % mean absolute price change at time k
SurvivorMAT = (1-lambda).*PDE;                                      % survivors at k 
end

MeanAbsPchange=MeanAbsPchange(2:end);                               % eliminating time 0 value
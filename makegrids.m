% Builds productivity and price grids in logs       

if finegrid==1                                      % use this for steady-state computation only
   nums = 101;                                      % number of cost grid points 
   NumPpoints = 101;                                % approximate number of price points; nump will be the actual
   numstddevs=5;                                    % number of standard deviations around midpoint
elseif finegrid==0                                  % use this for dynamics
   nums = 25;                                       % number of cost grid points 
   if accuracy>0
   NumPpoints = 31;                                 % number of price points on coarse grid
   numstddevs = 3;                                  % number of standard deviations around midpoint
   else
   NumPpoints = 31;                                 % number of price points on coarse grid
   numstddevs = 3;                                  % number of standard deviations around midpoint
   
   %NumPpoints = 25;                                 % number of price points on coarse grid       
   %numstddevs = 2.5;                                % number of standard deviations around midpoint
   end
elseif finegrid==-1                                 % use this for rep agent case
   nums = 1;                                        % number of cost grid points 
   NumPpoints = 3;                                  % number of price points on coarse grid
   numstddevs = 1;                                  % number of standard deviations around midpoint
end

SMAX=numstddevs*stdMC;                              % maximum cost
SMIN=-numstddevs*stdMC;                             % minimum cost
sstep=(SMAX-SMIN)/(nums-1);                         % distance between grid points
sgrid=SMIN:sstep:SMAX;                              % marginal cost grid (log)     

markup = log(epsilon/(epsilon-1));                  % optimal flexible price markup

if idioshocks == -1,                                % representative agent case
   nums=1;                                          % single productivity state
   SMAX=0;                                          % single grid point at zero log productivity
   SMIN=0;                                          % productivity=1 in levels
   sgrid=0;
   PMAX = markup + log(wflex) + gridSpread;         % need positive gridSpread for this to work       
   PMIN = markup + log(wflex) - gridSpread;          
else
%  Real price grid (LOGS)                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   PMAX = markup + log(wflex) + SMAX ;              % maximum price (based on flex price policy)
   PMIN = markup + log(wflex) + SMIN ;              % minimum price
   extraspread = gridSpread*(PMAX-PMIN);            % stretching out grid for rep agent mu=1 case: hit pstar=markup
   PMAX = markup + log(wflex) + SMAX + extraspread; % maximum price (based on flex price policy)
   PMIN = markup + log(wflex) + SMIN - extraspread; % minimum price
end

pstep = (PMAX-PMIN)/(NumPpoints-1);                 % NumPpoints price intervals
offset = log(mu)/pstep;                             % noninteger number of price steps caused by inflation
Pgrid=(PMIN:pstep:PMAX)';                           % price grid (log)
nump=length(Pgrid);                                 % number of price grid points
gridsize=nump*nums;                                 % total number of grid points 

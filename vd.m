clear STATEHISTORY; clear JumpHistory 
% Computes variance decomposition and Phillips curve regression

fprintf('\n')  
TT = 3*100;      % number of months 
INITCONDIT=0;

% SPECIFY MONEY SHOCK PROCESS for periods 1:TT
randn('state',0)
shocksize = jacstep;
scalefactor = 1/shocksize;
Rshocks = [shocksize *randn(1,TT)];  %#ok<NBRAK> %simulating random history
time1Rshock = Rshocks(1);

TFPshocks = zeros(1,TT);
time1TFPshock = TFPshocks(1);

distsim;
if phiPI > 0, compute_IRFs; else  compute_IRFsM; end    % run Taylor rule or money growth shock

if rem(TT,3)==0 
    
% convert to quarterly frequency
  C_pathQ = mean(reshape(C_path,3,TT/3));
  PI_pathQ = mean(reshape(PI_path,3,TT/3));
  d_R_pathQ = mean(reshape(d_R_path,3,TT/3));
  
  pchCQ = C_pathQ(2:end)./C_pathQ(1:end-1)-1;
  
% gdp growth and delfator inflation quarterly s.d. during 1984-2008
  data_std_cons_growth = 0.00510;  % okay:  cgg_mpr_qje.wf1
  data_std_infl        = 0.00246;  % okay:  cgg_mpr_qje.wf1 
  data_std_cons = 0.0090853;
  % computed as: std(dgdp), where dgdp = (gdp-hp_gdp)/hp_gdp and gdp is FRED II's series GDPC1: Real Gross Domestic Product, 1 Decimal
    

  scaleUp = data_std_infl/std(PI_pathQ);          % scale up to explain all observed inflation
  
  std_PI_pathQ  = std(PI_pathQ)*scaleUp;
  std_pchCQ     = std(pchCQ)*scaleUp;
  std_C_pathQ   = std(C_pathQ/Cbar)*scaleUp;

  VD_inflation   = std_PI_pathQ/data_std_infl  ;  % equals 1 by construction
  VD_cons_growth = std_pchCQ/data_std_cons_growth ;  
  VD_cons = std_C_pathQ/data_std_cons;
  
% 2SLS regression of consumption on inflation     
% first stage regression: inflation on exogenous shock 
  regressors = [ones(size(d_R_pathQ')) d_R_pathQ'];
  if rank(regressors)==size(regressors,2);
     B = regress(PI_pathQ',regressors);
  else
     error('Colinear regressors');
  end
  PI_projQ = B(1) + B(2)*regressors(:,2);
  
% second stage regression: output on predicted inflation     
  regressors = [ones(size(PI_projQ)) 4*log(PI_projQ)];
  if rank(regressors)==size(regressors,2);
     [B,BINT,R,RINT,STATS] = regress(log(C_pathQ'),regressors );
  else
     error('Colinear regressors');
  end

  
% Print output
  fprintf('\n')
  fprintf('100 x std dev of monetary shock                            : %0.4g \n',100*std(Rshocks)*scaleUp)
  fprintf('\n')
  fprintf('Model implied 100 x std of quarterly inflation             : %0.4g \n',100*std_PI_pathQ)
  fprintf('Actual 100 x std of quarterly deflator inflation 1984-2008 : %0.4g \n',100*data_std_infl)
  fprintf('Share of inflation variance due to monetary shocks         : %0.4g%% \n', 100*VD_inflation)
%   fprintf('\n')
%   fprintf('Model implied 100 x std of quarterly output growth         : %0.4g \n',100*std_pchCQ)
%   fprintf('Actual 100 x std of quarterly output growth 1984-2008      : %0.4g \n',100*data_std_cons_growth)
%   fprintf('Share of output variance due to monetary shocks            : %0.4g%% \n', 100*VD_cons_growth)
  fprintf('\n')
  fprintf('Model implied 100 x std of quarterly detrended output      : %0.4g \n',100*std_C_pathQ)
  fprintf('Actual 100 x std of quarterly detrended output 1984-2008   : %0.3g \n',100*data_std_cons)
  fprintf('Share of output variance due to monetary shocks            : %0.4g%% \n', 100*VD_cons)
  fprintf('\n')
  fprintf('Phillips curve regression: log(C_pathQ) = alpha + beta(4log(PI_projQ)) + eps \n')
  fprintf('\n')
  fprintf('Estimation method: 2SLS; Instrument for inflation: exogenous aggregate shock  \n')
  fprintf('Quarterly frequency (average of monthly simulated data) \n')
  fprintf('\n')
  fprintf('Slope coefficient beta                                     : %0.4g \n', B(2))
  fprintf('Standard error for slope coefficient                       : %0.4g \n', abs(B(2)-BINT(2,1))/2)
  fprintf('R2 of regression                                           : %0.4g ', STATS(1))
  fprintf('\n')
end

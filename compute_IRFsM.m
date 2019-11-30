% Computes dynamic paths as function of initial state and shocks (money rule)
% In particular computes impulse-response functions
% For simplicity, 'time' is same as vector index (i.e. starts at time=1)

% IN ORDER TO BE ABLE TO CONSIDER VERY SMALL SHOCKS, COMPUTE DEVIATIONS ONLY

% EXTRACTING VARIABLES DIRECTLY FROM KLEIN REPRESENTATION
% X = [Pdist, DP, mlag, V, C, PI, m]

% Extract exogenous shocks processes from STATEHISTORY
  d_R_path = STATEHISTORY(1,:);                                            % money growth process
  d_A_path = STATEHISTORY(2,:);                                            % TFP process

% Extract state dynamics from STATEHISTORY
  d_PhiHat_path = STATEHISTORY(nz+1:end-2,:);                              % history of LAGGED distributions
  PhiHat_path = Pdist(:)*ones(1,TT) + d_PhiHat_path;                       % by definition, PhiHat_t = Phi_t-1
  
  d_delta_wedge = STATEHISTORY(end-1,:);                                   % history of price dispersion
  delta_wedge = WeightPriceDispers + d_delta_wedge(2:end);                 % delta_wedge(t+1) is today's
  
  d_mlag_path = STATEHISTORY(end,:);                                       % history of lagged real money
  mlag_path = mbar + d_mlag_path ;                        
                                                                           
% Extract jump dynamics from JumpHistory
  d_V_path = JumpHistory(1:nV,:);                                          % value function history
  V_path = V(:)*ones(1,TT) + d_V_path;          

  d_C_path = JumpHistory(nV+1,:);                                          % consumption history
  C_path = Cbar + d_C_path;                                                

  d_PI_path = JumpHistory(nV+2,:);                                         % inflation history
  PI_path = mu + d_PI_path; 
  
  d_m_path = JumpHistory(nV+3,:);                                          % real money history
  m_path = mbar + d_m_path; 

% CONSTRUCT OTHER VARIABLES NOT INCLUDED IN KLEIN REPRESENTATION
% Nominal interest rate
  R_path = 1./(1 - nu*C_path.^gamma./m_path);
  R_path_a = R_path.^12 - 1;                                               % annualized net nominal interest rate 

% Real price level 
  p_path = 1./m_path;

% Real wage: from Labor supply equation (detrended by P)
  w_path  = chi*C_path.^gamma; 

% Labor hours path
  labor_path = C_path(1:end-1).*delta_wedge./exp(d_A_path(1:end-1));

% Nominal money growth path
  mu_path = PI_path.*m_path./mlag_path;

% Get optimal value function and optimal price
  Ms_path = NaN*ones(TT,nums); pStars_path = NaN*ones(TT,nums);     
  for time=1:TT
    Vt=reshape(V_path(:,time),nump,nums);
    [M_out pStar_out] = M_pStar(Vt, Pgrid, PI_path(time));
    Ms_path(time,:) = M_out;
    pStars_path(time,:) = pStar_out;
  end
  
% Adjustment gain and adjustment probability  
  D_path = kron(Ms_path',ones(nump,1)) - V_path;
  Lambda_path = adjustment(adjtype, V_path, D_path, ones(gridsize,1)*w_path, ksi, alpha, lbar, lam0);

  PhiTilde_path = NaN*zeros(gridsize,TT);         
  Rnowvec = NaN*zeros(nump^2,TT);         
  for time=1:TT
      Rnow = Rmatrix(nump, PI_path(time), pstep);
      PhiTildeNow = Rnow*reshape(PhiHat_path(:,time),nump,nums)*TRANSMAT';
      PhiTilde_path(:,time)= PhiTildeNow(:);       % PdistEroded history
      Rnowvec(:,time)=Rnow(:);
  end

% Now that we have constructed PhiTilde, we no longer need PhiHat.
  Phi_path = PhiHat_path;
  Phi_path(:,1) = [];  % here we are losing 1 period for Phi and everything computed from Phi
% So now we have history of Phi_t = distribution at time of production in period t.

% Compute mass of adjusters and fraction of adjusting firms
  Adjusters_path = Lambda_path.*PhiTilde_path;
  frac_adjusters_path = sum(Adjusters_path);

  Pchanges_path  = NaN*ones(nump,nums,TT); % needed for inflation decomposition later
  AvPchange_path = NaN*ones(1,TT);         % needed for inflation decomposition later
  for time=1:TT-1
    Pchanges_path(:,:,time) = ones(nump,1)*pStars_path(time,:) - Pgrid*ones(1,nums);
    AvPchange_path(time) = sum(sum(reshape(Pchanges_path(:,:,time),nump,nums)...
                           .*reshape(Adjusters_path(:,time),nump,nums)))./frac_adjusters_path(time);
  end
  
  ex_ante_real_interest_rate = R_path(1:end-1)-PI_path(2:end)-(Rss-mu);
  if INITCONDIT==0
    % the time 2 money shock is unexpected in period 1, so expected inflation equals steady-state
      ex_ante_real_interest_rate = [0 ex_ante_real_interest_rate(2:end)]; 
  end
  
% inflation decomposition
intensive_margin_path = zeros(1,TT-1); 
extensive_margin_path = zeros(1,TT-1);  
selection_effect_path = zeros(1,TT-1); 
dfracadj_path = zeros(1,TT-1); 
for time=1:TT
  infdecomp;
  intensive_margin_path(time) = intensive_margin; 
  extensive_margin_path(time) = extensive_margin;  
  selection_effect_path(time) = selection_effect; 
  dfracadj_path(time) = dfreqpchanges_dmu;
end



% Decomposes inflation into intensive, extensive, selection effects

%The steady state objects needed for this calculation are:
%   Pdisteroded: beginning-of-period distribution
%   lambda: the whole function lambda(p,a), nump x nums
%   Pdistans: matrix of desired LOG price adjustments (pstar-Pgrid), nump x nums 

format compact

%period in which money shock occurs - must be same as in "irf.m"
if exist('shocktime','var')
   if Rshocks(shocktime)>0
      d_shock = d_R_path(shocktime);
   elseif TFPshocks(shocktime)>0
      d_shock = d_A_path(shocktime);
   end
else
   d_shock = 1/scalefactor;
end
if d_shock == 0, d_shock = 1/scalefactor; end     % if there is no shock, take differences

% INFLATION IMPACT RESPONSE
dPI_dmu = (PI_path(time)-mu) / d_shock; 

% KLENOW-KRYVTSOV DECOMPOSITION
% inflation = freqpchanges*EPchange
dfreqpchanges_dmu = (frac_adjusters_path(time) - freqpchanges) / d_shock;
dAvPchange_dmu    = (AvPchange_path(time) - EPchange) / d_shock;

KKintensive = dAvPchange_dmu * freqpchanges;
KKextensive = dfreqpchanges_dmu * EPchange ;
KK_decomp = KKintensive + KKextensive ;

% NAKOV-COSTAIN DECOMPOSITION
% INFLATION = sum(sum(Pdistans.*(lambda - sum(sum(lambda .*Pdisteroded))).*Pdisteroded)) + ...
%    sum(sum(Pdistans.*Pdisteroded ))  *  sum(sum(lambda.*Pdisteroded )) ;
    
% intensive = sum(sum(Pdistans.*Pdisteroded )) ;                 % average DESIRED price change
% extensive = sum(sum(lambda.*Pdisteroded )) = freqpchanges  ;   % average frequency of adjustment
% selection = sum(sum(Pdistans.*(lambda - sum(sum(lambda .*Pdisteroded))).*Pdisteroded)) =
%           = extensive*(EPchange-intensive);                    % selection effect

d_Pchanges = (Pchanges_path(:,:,time) - Pdistans) / d_shock;
d_PhiTilde = (reshape(PhiTilde_path(:,time),nump,nums) - Pdisteroded) / d_shock;

d_intens = sum(sum(d_Pchanges.*Pdisteroded + Pdistans.*d_PhiTilde));  % element by element so derivative okay
d_extens = (frac_adjusters_path(time) - freqpchanges) / d_shock;

if adjtype==2
intensive_margin = dPI_dmu; 
extensive_margin = zeros(size(dPI_dmu)); 
else    
intensive_margin = d_intens*extensive; 
extensive_margin = d_extens*intensive; 
end
selection_effect = dPI_dmu - intensive_margin - extensive_margin;
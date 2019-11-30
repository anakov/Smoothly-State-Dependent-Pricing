function [lambda, POSSIBILITIES] = adjustment(adjtype, V, D, wage, ksi, alpha, lbar, lam0)
% Computes adjustment probability as a function of the gain from adjustment

L =  D./wage;                                          % Adjustment depends on D/w, value in units of labor time
switch adjtype
    case 0                                             % BASELINE SSDP MODEL 
      lambda = lbar./((1-lbar)*((alpha./L).^ksi)+lbar);    % probability of adjustment as a function of the state 
    case 2                                             % CALVO MODEL 
      lambda = lbar*ones(size(D));                     %   Constant probability of adjustment
    case 3                                             % GOLOSOV-LUCAS MODEL
      lambda = lamcontin(L,alpha);                         % probability of adjustment: 0 or 1
    case 4                                             % WOODFORD'S MODEL
      argexp=(alpha-L)*ksi;                                % argument in exponent function
      lambda = lbar./((1-lbar).*exp(argexp) + lbar);       % probability of adjustment as a function of the state 
    case 5                                             % DOTSEY-KING-WOLMAN MODEL
      lambda = lbar./((1-lbar)*((alpha./L).^ksi)+lbar);    % probability of adjustment as a function of the state 
end
if nargout==2
   if adjtype<3
     alpha = 0;                                        % No fixed cost is subtracted in Calvo and SSDP models
   elseif adjtype==5
     alpha = ExpectMenuCost(adjtype, D, wage, ksi, alpha, lbar, lam0); % Expected menu cost in DKW model 
   end
   POSSIBILITIES = V + (D - alpha*wage).*lambda;       % Continuation value in all models
end

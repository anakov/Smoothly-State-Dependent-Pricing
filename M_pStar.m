function [M, pStar] = M_pStar(V, Pgrid, PI)
% Computes value function maximum and argmax with quadratic interpolation on V

[nump,nums]=size(V);
M = NaN*ones(1,nums); 
pStar = NaN*ones(1,nums);   

[maxV, maxind] = max(V);                                    
if  min(maxind)==1 || max(maxind)==nump                                  % make sure solution is interior
    error('Corner solution. Increase gridSpread in Pidentity')
end

if  nums==1 && PI==1                                                     % if rep agent and zero inflation, then
    pStar = Pgrid(maxind);                                               % no interpolation: maximant is on the grid
    M  = maxV;                                                   
else
    localx = [Pgrid(maxind-1)';Pgrid(maxind)';Pgrid(maxind+1)'];         % 3 points on price grid around optimum
    XMATstack = [ones(size(localx));localx;localx.^2];                   % regressor matrix: const, grid and grid^2
    for col=1:nums
        XMAT = reshape(XMATstack(:,col),3,3);                            % regressor matrix at state "col"
        localV = V(maxind(col)-1:maxind(col)+1,col);                     % explained variable
        betacoeff = (XMAT'*XMAT)\XMAT'*localV;                           % regression coefficient
        pStar(col) = -0.5*betacoeff(2)/betacoeff(3);                     % maximant of interpolant
        M(col)  = [1 pStar(col) pStar(col)^2]*betacoeff;                 % maximum of interpolant
    end
end
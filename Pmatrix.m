function PmatOut = Pmatrix(pStar, Pgrid, pstep)

% Calculates matrix P which rounds stochastically to two grid points around the optimal price  

nump = length(Pgrid);
nums = length(pStar);
PmatOut = zeros(nump,nums);  % Pmat has zeros everywhere except around the diagonal, along the optimal price 

for col=1:nums                           
    OPTindHI=find(Pgrid>pStar(col),1);
    OPTindLO=max(OPTindHI-1,1);
    PmatOut(OPTindHI,col)= (pStar(col)-Pgrid(OPTindLO))/pstep;
    PmatOut(OPTindLO,col)= (Pgrid(OPTindHI)-pStar(col))/pstep;    
end                                                    


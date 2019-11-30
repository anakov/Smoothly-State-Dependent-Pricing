% Computes the stationary distributions of firms before and after shocks
% version 17 april 2008

PdistDIFF=inf;                                              % reset difference of Pdist
Pdistiter=0;                                                % reset counter
P = Pmatrix(pstar, Pgrid, pstep);                           % get distribution adjustment matrix

while (PdistDIFF>PdistDIFFtol && Pdistiter<1500)            % iterate to convergence of Pdist
    Pdistiter=Pdistiter+1;                                  % increment counter

    Pdisteroded = RECURSEMAT*Pdist*TRANSMAT';               % distribution after productivity shock and deflation
    PdistNew = (1-lambda).*Pdisteroded + P.*(ones(nump,nump)*(lambda.*Pdisteroded)); 
    
    PdistDIFF=gridsize*max(max(abs(PdistNew-Pdist)));       % sup norm normalized by gridsize
    Pdist=PdistNew;                                         % updating Pdist   
end

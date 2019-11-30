function R = taylor(Rss,phiR,Rlag,phiPI,PInow,mu,phiC,Cnow,Cbar,zRnow)
% Taylor rule
R = phiR*Rlag + (1-phiR)*( Rss + phiPI*(PInow-mu) + phiC*log(Cnow/Cbar) ) - zRnow;



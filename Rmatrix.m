function R = Rmatrix(nump, PI, pstep)

% Calculates the transition matrix R 
% (with the OFFSET DETERMINED by PI)

  R = sparse(nump,nump); 
  nowoffset = log(PI)/pstep;
  if nowoffset==0
      R = eye(nump);
  else
      remoffset = nowoffset - floor(nowoffset);
      
      R(1,1:ceil(nowoffset)) = 1;             
      
      startfirstdiag = [max([1 ; -floor(nowoffset)])  max([1 ; 1+ceil(nowoffset)])];         
      endfirstdiag = [min([nump-1 ; nump-ceil(nowoffset)])  min([nump ; nump+floor(nowoffset)])]; 
      
      startsecdiag = [max([2 ; 1-floor(nowoffset)])  max([1 ; 1+ceil(nowoffset)])];               
      endsecdiag = [min([nump ; nump-ceil(nowoffset)+1])  min([nump ; nump+floor(nowoffset)])];   
      
      R(startfirstdiag(1):endfirstdiag(1),startfirstdiag(2):endfirstdiag(2)) = ...
           remoffset*speye(nump-ceil(abs(nowoffset)));
       
      R(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) = ...
      R(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) + ...
           (1-remoffset)*speye(nump-ceil(abs(nowoffset)));
       
      R(nump,nump+floor(nowoffset)+1:nump) = 1; 
  end


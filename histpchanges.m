function masspchangesvec=histpchanges(Adjusters,logitprobMAT)
% Computes the histogram of price changes

[nump, nums] = size(logitprobMAT);

masspchanges=zeros(2*nump-1,nums);                                   % initialize mass in each bin with zero

for i = 1:nums                                                       % loop over states
adjprobMAT=Adjusters(:,i)*(logitprobMAT(:,i))';                      % distribution of adjusters at state i
  for j = 1-nump:1:nump-1
    masspchanges(j+nump,i)=sum(diag(adjprobMAT,j));                  % diagonal summing mass to price change bins
  end
end
masspchangesvec=sum(masspchanges,2);                                 % summing mass across states i

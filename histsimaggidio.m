% Simulates price histories with idiosyncratic productivity shocks 
% and aggregate monetary shocks

clear; load ssdp_tr; clc

lags = 60; 
training = 1200;      % need this high to reach ergodic distr. or prod. 
%months   = 245;      % number of months in Mackowiak et al 
months = 3*lags*30;   % number of months for shock identification

TT = training + months;
color='m';

FPS = 512;                   % firms per sector
S = 79;                    % number of sectors in Mackowiak et al 
N = FPS*S;                 % total number of firms           

rand('state',0)
randn('state',0)
rand1T  = randn(1,TT);
randTN  = rand(TT,N);
randNT  = rand(N,TT);
 
% scaling of shock to inflation volatility data
stdRshock = 0.00229;           % set to hit std of dlog(CPI) in Mackowiak

Rshocks = stdRshock*rand1T ;  % history of money shocks
time1Rshock = Rshocks(1);

% no aggregate TFP shock
TFPshocks = zeros(1,TT);          
time1TFPshock = TFPshocks(1);

scalefactor = 1; 
INITCONDIT = 0;                % starting from steady-state 

distsim;

if phiPI>0, compute_IRFs; else  compute_IRFsM; end

pricelevels = zeros(N,TT);
productivities = zeros(N,TT);
pricelevels(:,1)= 1;

fairy  = randNT;
sindex = round((nums+1)/2)*ones(N,1);   

CPI = cumprod(PI_path);
pstars_levels = exp(pStars_path') .* (ones(nums,1)*CPI);


for t=2:TT

if (rem(t,50)==1), clc, disp(sprintf('      Period       Total')),
    disp([t TT]), end  % report progress

% productivity shocks
[out, sindex] = max(cumsum(TRANSMAT(:,sindex))>=ones(nums,1)*randTN(t,:));

D = reshape(D_path(:,t),nump,nums);

eroded_current_real_prices = pricelevels(:,t-1)./CPI(t);

Dvec = interp2(1:nums, Pgrid, D, sindex', log(eroded_current_real_prices)); 

lambdavec=adjustment(adjtype, [], Dvec, w_path(t), ksi, alpha, lbar, lam0);

currentpstars = pstars_levels(:,t);

currentoptpricelevels = currentpstars(sindex);

pricelevels(:,t) = (fairy(:,t)>lambdavec).*pricelevels(:,t-1) + (fairy(:,t)<=lambdavec).* currentoptpricelevels ;

productivities(:,t)=exp(sgrid(sindex));

end

CPI  = CPI(:,training+1:end);
pricelevels = pricelevels(:,training+1:end);
productivities = productivities(:,training+1:end);
Rshocks = Rshocks(:,training+1:end);

sectorPriceIndex=NaN*zeros(S,length(CPI)); % sectoral price index
sectorProd=NaN*zeros(S,length(CPI)); % sectoral productivities

weights = [1./(1:FPS)/sum(1./(1:FPS))]'; % Zipf law for price weights 
s=0; 
for i=1:FPS:N 
s=s+1; 
sectorPriceIndex(s,:) = mean(pricelevels(i:i+FPS-1,:),1); 
sectorProd(s,:) = mean(productivities(i:i+FPS-1,:),1); 
secPriceZipf(s,:) = sum((pricelevels(i:i+FPS-1,:).* (weights*ones(1,months))),1); 
secProdZipf(s,:) = sum((productivities(i:i+FPS-1,:).* (weights*ones(1,months))),1); 
end 

CPIzipf =  mean(secPriceZipf);
%modgen = 100*[CPI/CPI(1) ; sectorPriceIndex/CPI(1) ]; 
modgen = 100*[CPIzipf/CPIzipf(1) ; secPriceZipf/CPIzipf(1)]; 

%agg = std(dlog(CPI)) ;
%sect = median(std(dlog(modgen(2:end,:)'))) ;
%sect_zipf = median(std(dlog(modgenZipf(2:end,:)'))) ;
%disp([100*std(dlog(CPI)) 0.1787])

%xlswrite('modgen_8zipf.xls',modgen)


%aggprod = log(mean(sectorProd))';
aggprod = log(mean(secProdZipf))';

for i=1:S
 clc;
 disp([i S]); 
 
% dlogsecprice = 100*dlog((sectorPriceIndex(i,:)))';
 dlogsecprice = 100*dlog((secPriceZipf(i,:)))';
  
%secprod = 100*log(sectorProd(i,:))';
 secprod = 100*log(secProdZipf(i,:))';
 prodinnov = secprod(2:end) - rho*secprod(1:end-1);
 
 aggprodinnov = aggprod(2:end) - rho*(aggprod(1:end-1));
 
 taylorshock = Rshocks(2:end)';
 
 
 X = [];    
 Y = dlogsecprice;
 x1 = prodinnov;
 x2 = taylorshock;
 %x2 = dlog(CPI)';
 x3 = aggprodinnov;
 
   for j=0:(lags-1)
      X = [X x1(lags-j:end-j,:) x2(lags-j:end-j,:) x3(lags-j:end-j)];
   end
 
Y = Y(lags:end,:);

%Y = standardize(Y); X = standardize(X);
X = [X ones(size(Y))];

B = REGRESS(Y,X);
Bs(:,i) = B(1:(end-1)); % excluding constant 

share(i) = sum(Bs(1:3:(3*lags),i).^2)*var(x1)/var(Y);
share2(i) = 1-sum(Bs(2:3:(3*lags),i).^2)*var(x2)/var(Y);
end

% figure(1)
% subplot(4,2,7)
% plot(cumsum(Bs(1:3:(3*lags),:)))
% title('Estimated responses to sectoral productivity shock')
% %ylim([0 1.5])
% subplot(4,2,8)
% plot(cumsum(Bs(2:3:(3*lags),:)))
% title('Estimated responses to aggregate interest rate shock')

mean(share2)


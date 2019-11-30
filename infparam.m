function mu = infparam(x)
% Loads a selection of money growth rates to choose from

switch x
    
case 'ACNiel'  
    mu = 0.999831;      % AC Nielsen inflation
case 'Domin'
    mu = 1.002822;      % Dominick's inflation
case 'GLlow'  
    mu = 1.0064^(1/3);  % GL07 baseline monthly money growth (2.58% annual)
case 'GLhigh' 
    mu = 1.21^(1/3);    % GL07 high monthly money growth (114.36% annual) 
case 'GanLow' 
    mu = 1.0037;        % Gagnon low (4.53% annual)
case 'GanMed' 
    mu = 1.0214;        % Gagnon medium (28.93% annual)
case 'GanHi'  
    mu = 1.0416;        % Gagnon high (63.08% annual)
case 'UShigh' 
    mu = 1.10^(1/12);   % US high in the 1970s (10% annual)
otherwise
    mu = str2num(x);      
end
    
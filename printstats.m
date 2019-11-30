% Prints out steady-state statistics
clc
disp(sprintf('\n'))
disp(sprintf('Frequency of price changes              : %0.1f%',100*freqpchanges))
disp(sprintf('Mean absolute price change              : %0.1f',100*MeanAbsPchange))
%disp(sprintf('Median absolute price change            : %0.1f',100*MedAbsPchange))
disp(sprintf('Std of price changes                    : %0.1f',100*STDpchange))  
disp(sprintf('Kurtosis of price changes               : %0.1f',KurtosisPchange))                                           
disp(sprintf('Fraction of price changes <5%%           : %0.1f%',100*fracSmallChanges))
if adjtype==3                                                        
   disp(sprintf('Menu cost as %% of revenue               : %0.1f%',McostinRev))
end
%disp(sprintf('Mean abs distance from p*               : %0.1f%',100*AvAbsDistFromPstar))
disp(sprintf('Mean loss in %% of frictionless profit   : %0.1f%',100*loss1))
disp(sprintf('Mean loss in %% of frictionless revenue  : %0.1f%',100*loss2))
disp(sprintf('Fit: Kolmogorov-Smirnov statistic       : %0.3f%',KSSTAT))
disp(sprintf('Fit: Euclidean distance                 : %0.3f%',EuclNorm))

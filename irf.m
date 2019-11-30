% Computes and plots impulse-response functions
clear STATEHISTORY; clear JumpHistory 
disp(sprintf('\n'))  
TT = 20;             
INITCONDIT=0;

% SET NOMINAL SHOCK IMPULSES
time1Rshock = 0;
time2Rshock = jacstep; 
time3Rshock = 0;

% SET TFP IMPULSES
time1TFPshock = 0;
time2TFPshock = 0;
time3TFPshock = 0;

% SPECIFY SHOCK IMPULSES for periods 1:TT
Rshocks   = [time1Rshock   time2Rshock   time3Rshock   zeros(1,TT-3)]; 
TFPshocks = [time1TFPshock time2TFPshock time3TFPshock zeros(1,TT-3)]; 

shocktime = 2;

if Rshocks(shocktime)>0
   scalefactor = abs(1/Rshocks(shocktime));    
elseif TFPshocks(shocktime)>0
   scalefactor = abs(1/TFPshocks(shocktime));  
end

distsim
if phiPI==0
   compute_IRFsM
elseif phiPI>0
   compute_IRFs
end
plot_IRFs




% Sets idiosyncratic shock and adjustment function parameters 
% minimum distance estimates on 31x25 grid (gridSpread 0.15)             
% check distance.m for criterion

switch adjtype
  case 0 % SDSP 0. SDSP 0. SDSP 0. SDSP 0. SDSP 0. SDSP 0. SDSP 0. SDSP 
        if version < 3 
                
        % baseline SSDP estimate
        lbar  = 0.110074147699003;
        alpha = 0.037210419907043;
        ksi   = 0.234597262202440;
        rho   = 0.900196054582284;
        stdMC = 0.155496755673145;
        
        % estimate on Dominicks
        % lbar  = 0.362500000000000;
        % alpha = 0.030407814930282;
        % ksi   = 0.040000000000000;
        % rho   = 0.900196054582284;
        % stdMC = 0.155496755673145;
        
        % estimate on 25x25 grid from JMCB paper
        % lbar  = 0.108928579454875;
        % alpha = 0.031083222345823;
        % ksi   = 0.293675662699625;    
        % rho   = 0.881156092163841;
        % stdMC = 0.147421080785644;
        end
  case 2 % CALVO 2. CALVO 2. CALVO 2. CALVO 2. CALVO 2. CALVO 2. CALVO 
        alpha = NaN;
        ksi   = NaN;
        lbar  =  0.10;
        %lbar  =  0.347; % Dominicks
        if version==1     
        rho   = 0.854023513626164;
        stdMC = 0.163449540209002;
 
        % estimate on 25x25 grid from JMCB paper 
        % rho   = 0.857625081201242;
        % stdMC = 0.165363367309271;
        end
  case 3 % FIXED MENU COST 3. FIXED MENU COST 3. FIXED MENU COST 
        lbar  = NaN;
        ksi   = NaN;
        if version==1
        % main estimate for JME paper    
        alpha = 0.066520066968469;
        rho   = 0.827959555912745;
        stdMC = 0.137529687188112;

        % estimate on Dominick's
        % alpha = 0.0063;
        % rho   = 0.900196054582284;
        % stdMC = 0.155496755673145;
        
        % estimate on 25x25 grid from JMCB paper
        % alpha = 0.0630548117601; 
        % rho   = 0.8468875723353; 
        % stdMC = 0.1445805243468;
        end
  case 4 % WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 4. WOODFORD 
        if version==1        
        lbar  = 0.094511020919308;
        alpha = 0.061077348715626;
        ksi   = 1.333516771762109;
        rho   = 0.857490512850195;
        stdMC = 0.179512443393469;
        % estimate on 25x25 grid from JMCB paper
        % lbar  = 0.0946;
        % alpha = 0.0609;
        % ksi   = 1.3341;
        % rho   = 0.8596;
        % stdMC = 0.1805;
        end
  case 5 % SMC 5. SMC 5. SMC 5. SMC 5. SMC 5. SMC 5. SMC  
        if version==1                  
        lbar  = 0.110015409129278;
        alpha = 0.037313524555662;
        ksi   = 0.235115626924701;
        rho   = 0.900293327913178;
        stdMC = 0.155388363572221;  
        % estimate on 25x25 grid from JMCB paper
        % lbar  = 0.108579977283206;
        % alpha = 0.034863241298866;
        % ksi   = 0.231778803383148;
        % rho   = 0.918308821866637;
        % stdMC = 0.160007325336373;
        end
  case 8 % SSDP L(0)>0 8. SSDP L(0)>0 8. SSDP L(0)>0 8. SSDP L(0)>0 
        if version==1                              
        lbar  = 0.108992820185848;
        alpha = 0.040606382073371;
        ksi   = 0.226623569306285;
        rho   = 0.903399950390407;
        stdMC = 0.155676640466923;
        lam0  = 0.001107331289855;
        end
end
if adjtype~=8, lam0=NaN; end



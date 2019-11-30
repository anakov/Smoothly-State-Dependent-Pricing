% Reports convergence statistics

function Y = convergence_report(adjtype,finegrid,fzeroiter,VDIFF,Viter) 
clc
disp(sprintf('COMPUTING STEADY-STATE:',adjtype))
disp(sprintf('\n'))  
disp(sprintf('Model      : %d',adjtype))
disp(sprintf('Grid type  : %d',finegrid))
disp(sprintf('\n'))  
disp(sprintf('Fzero iter : %d',fzeroiter))
disp(sprintf('V DIFF     : %0.3d',VDIFF))
disp(sprintf('V iter     : %d',Viter))
Y=[];
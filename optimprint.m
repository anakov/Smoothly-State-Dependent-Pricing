function stop = optimprint(x,optimValues,state)
% Prints current parameter vector at each iteration of the estimation procedure
%
%   STOP = OPTIMPRINT(X,OPTIMVALUES,STATE) prints the current point, X
%
%   Example:
%   Create an options structure that will use OPTIMPRINT
%   as the print function
%       options = optimset('OutputFcn',@optimprint);
%
%   Pass the options into an optimization problem to view the plot
%       fminbnd(@sin,3,10,options)

stop = false;
switch state
case 'iter'
      x = x(:);
      disp(sprintf('Parameter vector:'))  
      disp(x)
      disp(sprintf('Distance metric:'))  
      disp(optimValues.fval)
      if exist('parammat.mat','file'), load parammat, else paramvec = []; end
      x = [x ; optimValues.fval];
      paramvec = [paramvec x];
      save parammat paramvec 
end
end


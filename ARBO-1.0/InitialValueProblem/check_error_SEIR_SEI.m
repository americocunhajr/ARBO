% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
% Initial value problem: check_error_SEIR_SEI.m
%
% This is a main file for check the rhs_SEIR_SEI.m 
% function constraints
%
%
% Inputs:
%   t: time                                 - double
%   y: state vector                         - double array (5x1)
%   param: parameters vector                - double array (?x1)
%   rhs_SEIR_SEI: SEIR-SEI equations file   - .m function file
%
% Output:
%   log                                     - interface messages
% -----------------------------------------------------------------
% programmer: Michel Tosin
%             michel.tosin@uerj.br
%
% number of lines: 330
% last update: Jun 4, 2021
% -----------------------------------------------------------------

clc
clear
close all

% display program header
% -----------------------------------------------------------------
disp('                                                            ')
disp('============================================================')
disp('             ARBO - Arbovirus Modeling and                  ')
disp('          Uncertainty Quantification Toolbox                ')
disp('                                                            ')
disp('            by Eber Dantas, Michel Tosin,                   ')
disp('       Americo Cunha Jr and Rebecca E. Morrison             ')
disp('                                                            ')
disp('     This is a package for simulation and analysis          ')
disp('            of arbovirus nonlinear dynamics.                ')
disp('============================================================')
disp('                                                            ')
disp(' -----------------------------------------------------------')
disp(' ++++++++++++++++++ rhs_SEIR_SEI checking +++++++++++++++++ ')
disp(' -----------------------------------------------------------')
% -----------------------------------------------------------------

% Case 1: length(param) < 7 
% -----------------------------------------------------------------
t = 1;                         % time
y = [999 0 1 0 1 0 0 1];       % state vector
param = [1000 0.5*ones(1,5)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 1: length(param) < 7')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 2: length(param) > 7
% -----------------------------------------------------------------
param = [1000 0.5*ones(1,7)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 2: length(param) > 7')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 3: mod(param(1),1) ~= 0
% -----------------------------------------------------------------
param = [1000.25 0.5*ones(1,6)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 3: mod(param(1),1) ~= 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 4: param(1) < 0
% -----------------------------------------------------------------
param = [-1000 0.5*ones(1,6)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 4: param(1) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 5: param(2) < 0
% -----------------------------------------------------------------
param = [1000 -0.5 0.5*ones(1,5)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 5: param(2) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 6: param(3) < 0
% -----------------------------------------------------------------
param = [1000 0.5 -0.5 0.5*ones(1,4)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 6: param(3) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 7: param(4) < 0
% -----------------------------------------------------------------
param = [1000 0.5 0.5 -0.5 0.5*ones(1,3)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 7: param(4) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 8: param(5) < 0
% -----------------------------------------------------------------
param = [1000 0.5*ones(1,3) -0.5 0.5 0.5];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 8: param(5) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 9: param(6) < 0
% -----------------------------------------------------------------
param = [1000 0.5*ones(1,4) -0.5 0.5];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 9: param(6) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 10: param(7) < 0
% -----------------------------------------------------------------
param = [1000 0.5*ones(1,5) -0.5];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 10: param(7) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 11: length(y) < 8
% -----------------------------------------------------------------
y = [999 0 1 0 1 0 0];         % state vector
param = [1000 0.5*ones(1,6)];  % parameters vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 11: length(y) < 8')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 12: length(y) > 8
% -----------------------------------------------------------------
y = [999 0 1 0 1 0 0 1 1];     % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 12: length(y) > 8')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 13: y(1) < 0
% -----------------------------------------------------------------
y = [-999 0 1 0 1 0 0 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 13: y(1) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 14: y(2) < 0
% -----------------------------------------------------------------
y = [998 -1 1 0 1 0 0 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 14: y(2) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 15: y(3) < 0
% -----------------------------------------------------------------
y = [999 0 -1 0 1 0 0 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 15: y(3) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 16: y(4) < 0
% -----------------------------------------------------------------
y = [999 0 1 -1 1 0 0 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 16: y(4) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 17: y(5) < 0
% -----------------------------------------------------------------
y = [999 0 1 0 -1 0 0 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 17: y(5) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 18: y(6) < 0
% -----------------------------------------------------------------
y = [999 0 1 0 1 -0.1 0 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 18: y(6) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 19: y(7) < 0
% -----------------------------------------------------------------
y = [999 0 1 0 1 0 -0.1 1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 19: y(7) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 20: y(8) < 0
% -----------------------------------------------------------------
y = [999 0 1 0 1 0 0 -1];      % state vector

try  
    rhs_SEIR_SEI(t,y,param)
catch ME
    disp(' ')
    disp('Case 20: y(8) < 0')
    disp(['t: ',num2str(t)])
    disp(['y: ',num2str(y)])
    disp(['param: ',num2str(param)])
    disp(' ')
    disp(ME.message)
    disp(' ')
    disp('Error log checked!')
    disp(' ')
end
% -----------------------------------------------------------------

% Case 21: valid input
% -----------------------------------------------------------------
y = [999 0 1 0 1 0 0 1];      % state vector
 
dy = rhs_SEIR_SEI(t,y,param);

disp(' ')
disp('Case 21: valid input')
disp(['t :',num2str(t)])
disp(['y: ',num2str(y)])
disp(['param: ',num2str(param)])
disp(' ')
disp(['dy: ',num2str(dy')])
disp(' ')
disp('Valid input checked!')
disp(' ')
% -----------------------------------------------------------------

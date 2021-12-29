% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
% Calibration: ObjFunIC_SEIR_SEI.m
%
% This file is used especifically for the Trust Region algorithm.
% This file supplies the F(p,xdata) for the ydata comparisson.
% F(p,xdata) is extracted from the ode45 responde vector.
% Check main Trust Region associated file for futher explanation.
%
% The quantity of interest here is the number of new human cases 
% in a SEIR-SEI model. For more details about the model, check 
% the routines in the module Initial Value Problem. 
%
% Inputs:
%   param: parameters vector                   - double array (7x1)
%   IC: initial conditions vector              - double array (7X1)
%   tspan: time interval                       - double array (?x1)
%   tu: time unit reference                    - double
%   time: the restrict interval in tu units    - double array (?x1)
%   rhs_SEIr_SEI: SEIR-SEI equations file  - .m function file
%
% Output:
%   F: number of new cases in time and IC sums - double array (?x1)
% -----------------------------------------------------------------
% programmers: Eber Dantas
%              eberdantas@coppe.ufrj.br
%
%              Michel Tosin
%              michel.tosin@uerj.br
%
%              Americo Cunha
%              americo.cunha@uerj.br
%
% number of lines: 9
% last update: Jun 9, 2021
% -----------------------------------------------------------------

% Function
% -----------------------------------------------------------------
function F = ObjFunIC_SEIR_SEI(time,param,IC,tspan,tu)

%% ODE solver Runge-Kutta45
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
[~, y] = ode45(@(t,y)rhs_SEIr_SEI(t,y,param),tspan,IC,opt);

% time series of cumulative number of infected humans (number of ppl)
C  = y(:,7);

% discrete number of new cases
dt = tspan(2)-tspan(1); % stpan time step
h = round(tu/dt);       % step of convertion
Ctu = C(1:h:end);       % extract the C in the tu instants

% Number of new infectious humans
NewCases = [Ctu(1); ( Ctu(2:end) - Ctu(1:end-1) )]; 

% human initial conditions sum
%SumH = SH0+EH0+IH0 (= N - R0)
SumH = y(1,1) + y(1,2) + y(1,3);  

% vectors initial ocnditions sum
% SumV = SV0+EV0+IV0 (= 1)
SumV = y(1,4) + y(1,5) + y(1,6);  

% Quantity of interest in time instants and Sums   
F = [NewCases(time(1:end-2)); SumH; SumV];

end
% -----------------------------------------------------------------

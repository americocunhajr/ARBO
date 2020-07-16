%% Function output for the TRR (IC) algorithm on SEIR-SEI ZIKAV system 
% -----------------------------------------------------------------
% This file is used especifically for the Trust Region algorithm (with IC).
% This file supplies the F(p,xdata) for the ydata comparisson.
% F(p,xdata) is extracted from the ode45 responde vector.
% Check main Trust Region associated file for futher explanation.

%   day               = the day corresponding to each empirical data
%   case_name         = case configuration. Need to be passed down to the auxiliary function 
%   variable_TRRparam = TRRparameters varied by the algorithm (includes system parameters and/or initial conditions
%   fixed_ TRRparam   = TRRparameters fixed
%   IC                = vector of initial conditions
%   tspan             = time interval. Used on the ode's integration
%   dt                = time step of the ode's integration
%   Ndt               = number of time steps. Corresponds to the lenght of the time vector
%
% ----------------------------------------------------------------- 
%  programmers: Eber Dantas
%               eber.paiva@uerj.br
%
%               Michel Tosin
%               michel.tosin@uerj.br
%
%               Americo Cunha
%               americo.cunha@uerj.br
%
%  last update: May 2, 2017
% -----------------------------------------------------------------

%% Function
% -----------------------------------------------------------------
function F = TRR_FunctionOutput_SEIR_SEI(day,case_name,variable_TRRparam,fixed_TRRparam,tspan,dt,IC)

%% Adjusting IC vector

% Default TRRparam vector: [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0 C0]
% IC = [SH0 EH0 IH0 SV0 EV0 IV0 C0]
    % RH0 is omitted to optimize the procedure

if strcmp(case_name,'case_config1');

    % variable_TRRparam = [IV0];
    IC(6) = variable_TRRparam(1);
    
elseif strcmp(case_name,'case_config2');

    % variable_TRRparam = [EV0 IV0];
    IC(5) = variable_TRRparam(1);
    IC(6) = variable_TRRparam(2);

elseif strcmp(case_name,'case_config3');

    % variable_param = [bH bV EV0 IV0]
    IC(5) = variable_TRRparam(3);
    IC(6) = variable_TRRparam(4);

elseif strcmp(case_name,'case_config4');

    % variable_TRRparam = [bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0]
    IC    = [variable_TRRparam(7:12) IC(7)];

elseif strcmp(case_name,'case_config5');

    % variable_TRRparam = [bH aH yH bV aV dV EH0 IH0  EV0 IV0]
    SH0 = fixed_TRRparam(1) - variable_TRRparam(7) - variable_TRRparam(8) - 29639;
    SV0 = 1 - variable_TRRparam(9) - variable_TRRparam(10);
    IC    = [SH0 variable_TRRparam(7:8) SV0 variable_TRRparam(9:10) IC(7)];
    
end

%% ODE solver Runge-Kutta45
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
[time y] = ode45(@(t,y)TRR_rhs_SEIR_SEI(t,y,case_name,variable_TRRparam,fixed_TRRparam),tspan,IC,opt);

% y = [SH EH IH SV EV IV C] is the solution vector
    % RH is omitted to optimize the procedure

% time series of cumulative number of infected humans (number of ppl)
C  = y(:,7);

% SumH = SH0+EH0+IH0 (= N - R0)
%SumH = y(1,1) + y(1,2) + y(1,3);  

% SumV = SV0+EV0+IV0 (= 1)
%SumV = y(1,4) + y(1,5) + y(1,6);  


% discrete number of new cases
NewCases = C(1);
for i = 2:52
    NewCases(i) = C( (7*i-7)/dt +1 ) - C( (7*(i-1)-7)/dt +1 );
end

%% Position on the time vector corresponding to the values of 'day'

% Comparing C
    % Subtract 'day' value by 'time(1)' then divide by the timestep to find the
    % number of steps on the time matrix between time(1) and the day value. Add
    % 1 to count the first element.
%Tposition = (day - time(1))/dt + 1; 
       
% Comparing NewCases
    % The NewCases vector has the same number of elements than the day vector
    % Subtract two elements, they correspond to the Sums comparisson, not the NewCases one.
Tposition = day(1:end)/7; 


%% F(p,xdata). Output value

% Comparing NewCases and Sums
F = [NewCases(Tposition)]; %SumH SumV];

end
% -----------------------------------------------------------------

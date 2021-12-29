% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
% Calibration: TRR_main_SEIR_SEI_example2.m
%
% This is the main file for the Trust-Region-Reflective lsqcurvefit 
% for the SEIR-SEI ZIKAV Brazil outbreak. More ifnormation about 
% the model can be found in the initial problem value module files.
%
%   Procedure:
%   The function lsqcurvefit(fun,inital_guess_param,xdata,ydata) provides
%   values for the varying parameters that best approximate F(p,xdata) to
%   the ydata points on a least-squared-differences sense
%   - fun is the function that supply F(p,xdata)
%   - p is the vector containing the varying parameters.
%   - xdata is the data that correspond to the domain values
%   - F(p,xdata) is the function values that must approximate to ydata
%   - ydata is the fel data
%   - initial_guess_param is the initial guess of p that must be supplied
%   to the algorithm
%
%   When using a ode45 function, a second function must be made that provides
%   F(p,xdata). F(p,xdata) will not be the full result of ode45
%   integration. This second function must evaluate the 'y' vector response 
%   of ode45 and pick the F(p,xdata) values on the 'y' vector that 
%   correspond to the ydata values for comparison.
%   A third function is required to integrate through ode45, the
%   system defining function of the integration.
%   The associated fixed and varying parameters must be specified in each
%   function. 
%
% This code uses TRR_rhs_SEIR_SEI.m to define the ODE system and 
% the TRR_QoI_SEIR_SEI.m file to calculate the model output.
%
% In this example is calibrated the follow model parameters:
%
%   betaH  = vector-to-human transmission rate   (days^-1)
%   alphaH = human latent rate                   (days^-1)
%   gamma  = human recovery rate                 (days^-1)
%
%   betaV  = human-to-vector transmission rate   (ind^-1days^-1)
%   alphaV = vector latent rate                  (days^-1)
%   delta  = vector birth/mortality rate         (days^-1)
%
% and the follow initial conditions sums:
%
%   SumH   = sum of the human SEI initial conditions (individuals)
%         % RH0 is omitted to optimize the procedure 
%   SumV   = sum of the human SEI initial conditions (dimensionless)
%
% The conditions simulated by this file was estrated 
% from the Brazil's zika outbreak discrebed in 
% E. Dantas, M. Tosin and A. Cunha Jr,
% Uncertainty quantification in the
% nonlinear dynamics of Zika virus, 
% v3, 2020.
%
% Inputs:
%   param: parameters vector                    - double array (7x1)
%   IC: initial conditions vector               - double array (7X1)
%   tspan: time interval                        - double array (?x1)
%   rhs_SEIr_SEI: model equations file          - .m function file
%   ObjFunIC_SEIR_SEI: model output file        - .m function file
%
% Outputs:
%   TRRparam: calibrated set of parameters      - double array (6x1)
%   R_nought: basic reproduction number         - double
%   figure 1: model state in time               - inplace figure
%   figure 2: number of new cases in time       - inplace figure
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
% number of lines: 114
% last update: Jun 9, 2021
% -----------------------------------------------------------------


clc
clear
close all


% Data - Brazil
% -----------------------------------------------------------------
% Load 2016 data
ProbData = load('2016zika_Prob250117.dat');

% Data time vector
timeData = ProbData(1,:)';

% New cases data
NewCasesData = ProbData(2,:)';

% Cumulative cases data
CumCasesData = cumsum(NewCasesData);

% Load 2015 data - Cumulative number of cases
C2015 = 29639;      % until EW 49
% -----------------------------------------------------------


% parameters and initial conditions [USER INPUT]
% -----------------------------------------------------------      
% human population size (number of individuals)
N = 206e6; % Brazil population in 2016          

% transmission rates
betaH = 1/11.3; % vector-to-human (days^-1) 
betaV = 1/8.6;  % human-to-vector (ind^-1days^-1) 

% latent periods (days)
TalphaH = 5.9;
TalphaV = 9.1;
    
% latent rates (days^-1)
alphaH =  1/TalphaH;
alphaV =  1/TalphaV;
        
% human recovery period (days)
Tgamma = 7.9;

% human recovery rate (days^-1)
gamma  = 1/Tgamma;
    
% Mosquito lifespan (days)
Tdelta = 11;
    
% Mosquito birth/mortality rate (days^-1)    
delta = 1/Tdelta;

% initial conditions    
C0 = CumCasesData(1); % initial cumulative infectious humans (number of individuals)

IH0 = C0;                  % initial infectious humans  (number of individuals)
EH0 = IH0;                 % initial exposed humans     (number of individuals)
RH0 = C2015;               % initial recovered humans   (number of individuals)
SH0 = N - EH0 - IH0 - RH0; % initial susceptible humans (number of individuals)

IV0 = 2.2e-4;              % initial infectious mosquitoes      (dimensionless)
EV0 = IV0;                 % initial exposed mosquitoes         (dimensionless)
SV0 = 1 - EV0 - IV0;       % initial susceptible mosquitoes     (dimensionless)

SumHData = N - RH0;
SumVData = 1;
     
% -----------------------------------------------------------

% TRR procedure
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- Running TRR Algorithm --- ');
disp(' ');
disp(' ');


% time interval of analysis
   t0 = 7;                  % initial time (days)
   t1 = 364;                % final time   (days)
   dt = 1  ;                % time steps   (days)
tspan = t0:dt:t1;           % interval of analysis
   tu = 7;

% Output function
    % This function must supply the F(p,xdata) for the ydata comparisson

F = @(Param,time)ObjFunIC_SEIR_SEI(time,[N Param(1:6)],[Param(7:end) C0],tspan,tu); 


% TRR Algorithm

    % Choose bounds
    % Choose time series comparisson. Remember to change on the TRR_FunctionOutput accordingly.
    % TRRparam display the new parameters in the order of TRRparam0
       
paramlb = [1/16.3 1/12 1/8.8 1/11.6 1/10 1/21];      % lower bound result of variable parameters
paramub = [ 1/8   1/3   1/3  1/6.2  1/5  1/11];      % upper bound result of variable parameters
TRRparam0 = [betaH alphaH gamma betaV alphaV delta]; % initial guess for the parameters

IClb = [0.9*N  0   0   0.99 0 0];                    % lower bound result of variable IC
ICub = [  N   1e4 1e4 0.999 1 1];                    % upper bound result of variable IC
TRRIC0 = [SH0 EH0 IH0 SV0 EV0 IV0];                  % initial guess for the IC

lb = [paramlb IClb]; % lower bound for the TRR variable inputs
ub = [paramub ICub]; % upper bound for the TRR variable inputs

TRRguess = [TRRparam0 TRRIC0];
    
options = optimset('Display','iter','Algorithm','trust-region-reflective','MaxIter',1000,'TolFun',1e-10);
%    options = optimoptions('lsqcurvefit','Display','iter','Algorithm','trust-region-reflective','MaxIterations',100000,'MaxFunctionEvaluations',%500000,'FunctionTolerance',1e-007,'StepTolerance',5e-007,'FiniteDifferenceType','central');
            % Display: show process/result. 'iter' whill show the steps.
            % Algorithm: choose algorithm. LM is default. LM can't use lb and up
            % MaxIterations: maximum number of iterations. Default is 400
            % MaxFunEvals: maximum number of function evaluations allowed. Default is 100*NumberOfVariables.
    
            
% Comparing NewCases
TRRresult = lsqcurvefit(F,TRRguess,[timeData; 400; 500],[NewCasesData; SumHData; SumVData],lb,ub,options);

% -----------------------------------------------------------


% Integrate the new initial value problem
% -----------------------------------------------------------
disp(' ');
disp(' --- integration of the nonlinear dynamics with new parameters and IC --- ');
disp(' ');
disp('    ... ');
disp(' ');

TRRNewCases = F(TRRresult,[timeData; 400; 500]);
TRRCumCases = cumsum(TRRNewCases);

% basic reproduction number
R0 = ( TRRresult(1)*TRRresult(4)*TRRresult(5) )/( TRRresult(3)*(TRRresult(5)+TRRresult(6))*TRRresult(6) );


% plot Cumulative cases
figure(1)
fig1(1) = scatter(timeData,TRRCumCases(1:end-2));
hold on
fig1(2) = scatter(timeData,CumCasesData);
hold off

    % plot labels
     title('Cumulative human cases (SEIR-SEI)');
    xlabel('time (weeks)'                 );
    ylabel('number of individuals'        );

    % set plot settings
    set(gca,'FontSize',18);
    set(fig1(1),{'MarkerEdgeColor'},{'b'});
    set(fig1(2),{'MarkerEdgeColor'},{'m'});
    set(fig1(1),{'LineWidth'},{3});
    set(fig1(2),{'LineWidth'},{3});

    % legend
    leg = {'model';'data'};
    legend([fig1(1),fig1(2)],leg,'FontSize',15,'location','northeast');
    
    saveas(figure(1),'fig_TRR_example2_C.png')
    
% plot NewCases (per week) of SEIR-SEI model
figure(2)
fig2(1) = stem(timeData,TRRNewCases(1:end-2));
hold on
fig2(2) = stem(timeData,NewCasesData);
hold off

    % plot labels
     title('New Cases per week (SEIR-SEI)');
    xlabel('time (weeks)'                 );
    ylabel('number of individuals'        );

    % set plot settings
    set(gca,'FontSize',18);
    set(fig2,{'Color'},{'b';'m'});
    set(fig2,{'LineWidth'},{3;3});

    % legend
    leg = {'model';'data'};
    legend(fig2,leg,'FontSize',15,'location','northeast');
    
    saveas(figure(2),'fig_TRR_example2_NC.png')
% -----------------------------------------------------------------   

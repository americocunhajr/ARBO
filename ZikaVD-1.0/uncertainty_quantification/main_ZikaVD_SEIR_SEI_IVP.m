
%Main file for the SEIR-SEI ZikaV Brazil outbreak
% -----------------------------------------------------------
% Calculations are made on a day time scale and plotting on a week time
% scale
%
%   N = human population size (number of ppl)
%
%   SH = susceptible humans (number of ppl at t)
%   EH = humans incubating (number of ppl at t)
%   IH = infectious humans (number of ppl at t)
%   RH = recovered humans (number of ppl at t)
%
%   bH = vector-to-human transmission rate (days^-1)
%   aH = incubation human frequency (days^-1)
%   yH = infectious human frequency (days^-1)
%
%   C  = cumulative number of infectious humans (number of ppl at t)
%
%
%   SV = proportion of susceptible vector (adimensional)
%   EV = proportion of vectors incubating (adimensional)
%   IV = proportion of infectious vectors (adimensional)
%
%   bV = human-to-vector transmission rate (days^-1)
%   aV = incubation vector frequency (days^-1)
%   dV = vector lifespan frequency (days^-1)
% 
% ----------------------------------------------------------------- 
%  programmer: Eber Dantas
%              eber.paiva@uerj.br
%              
%              Michel Tosin
%              michel.tosin@uerj.br
%
%              Americo Cunha
%              americo@ime.uerj.br
%
%  last update: March 1, 2019
% -----------------------------------------------------------------


% Header
% -----------------------------------------------------------
clc
clear all
close all

disp(' '); 
disp(' --- Zika SEIR-SEI model simulation for Brazil outbreak --- ');
disp(' ');
disp('    ... ');
disp(' ');
% -----------------------------------------------------------


% Data - Brazil
% -----------------------------------------------------------------
disp(' '); 
disp(' --- loading data --- ');
disp(' ');
disp('    ... ');
disp(' ');

% Load 2016 data
ProbData = load('2016zika_Prob250117.dat');

ConfData = load('2016zika_Conf250117.dat');

% Data time vector
day = ProbData(1,:)*7;
%day = ConfData(1,:)*7;

% New cases data
NewCasesData = ProbData(2,:);
%NewCasesData = ConfData(2,:);

% Cumulative cases data
CData = cumsum(NewCasesData);

% Load 2015 data - Cumulative number of cases
%C2015 = 27127;     % until EW 45
C2015 = 29639;      % until EW 49

%C2015 = 443502;     % inferior PAHO estimation     
%C2015 = 1301140;    % superior PAHO estimation
% -----------------------------------------------------------


% Physical parameters and initial conditions
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters and initial conditions --- ');
disp(' ');
disp('    ... ');
disp(' ');
       
        
% Human population size (number of ppl)
    N = 206e6;          


% Intrinsic incubation period (days)
    TaH = 12;
    %TaH = 5.9-1.5;
        
% Extrinsic incubation period (days)
    %TaV = 15;                
    TaV = 10;
        
% Human infectious period (days)
    TyH = 3;
    %TyH = 9.9-TaH-1.5;
        
% Mosquito lifespan (days)
    TdV = 21;        
    %TdV = 7;
    %TdV = 14;
    %TdV = 30;
    
    
% Rates and inverses
    aH  = 1./TaH;                   % intrinsic incubation ratio (days^-1)

    aV  = 1./TaV;                   % extrinsic incubation ratio (days^-1)

    yH  = 1./TyH;                   % human infectious ratio (days^-1)
       
    dV  = 1./TdV;                   % mosquito death rate (lifespan inverse) (days^-1)

    
% Transmission rates
    bH  = 1/10.40;               % vector-to-human transmission rate (days^-1)
    %bH  = 0.58;
    bV  = 1/7.77 ;               % human-to-vector transmission rate (days^-1)
    %bV  = 0.59;

% Initial conditions
    %IH0 = NewCasesData(1);          % initial number of infectious humans
    IH0 = 10000;    
    %EH0 = IH0;                      % initial number of incubating humans
    EH0 = 6827;    
    %SH0 = N - IH0 - EH0 - C2015;    % initial number of susceptible humans
    SH0 = 205953534;    
    RH0 = C2015;                        % initial number of recovered humans
        
    %C0  = IH0;                      % initial cumulative number of infectious humans
    C0 = NewCasesData(1);
        
    %IV0 = 0.00022;                  % initial proportion of infectious vectors
     IV0 = 0;  
    EV0 = 0.000414;                      % initial proportion of incubating vectors
        
    %SV0 = 1 - IV0 - EV0;            % initial proportion of susceptible vectors
    SV0 = 0.999; 
% -----------------------------------------------------------


% Integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nonlinear dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');


% initial time of analysis
t0 = 7;

% final time of analysis
t1 = 365;                   % analysis starting at january 2016

% value of time steps
dt = 1;

% interval of analysis
tspan = t0:dt:t1;
Ndt   = length(tspan);      % number of time steps

% physical parameters vector
phys_param = [N bH aH yH bV aV dV];

% Initial conditions vector
IC = [SH0 EH0 IH0 RH0 SV0 EV0 IV0 C0];

% ODE solver Runge-Kutta45
        % ODE solver optional parameters
            %opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
[time y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,phys_param),tspan,IC);


% time series of susceptible humans (number of ppl)
SH = y(:,1);

% time series of incubating humans (number of ppl)
EH = y(:,2);

% time series of infectious humans (number of ppl)
IH = y(:,3);

% time series of recovered humans (number of ppl)
RH = y(:,4);


% time series of proportion of susceptible vectors (number of vectors)
SV = y(:,5);

% time series of proportion of incubating vectors (number of vectors)
EV = y(:,6);

% time series of proportion of infectious vectors (number of vectors)
IV = y(:,7);


% time series of cumulative number of infected humans (number of ppl)
C  = y(:,8);

% discrete number of new cases
NewCases = C(1);
for i = 2:52
    NewCases(i) = C( (7*i-7)/dt +1 ) - C( (7*(i-1)-7)/dt +1 );
end

toc
% -----------------------------------------------------------


% -----------------------------------------------------------
post_ZikaVD_SEIR_SEI_IVP
% -----------------------------------------------------------


    

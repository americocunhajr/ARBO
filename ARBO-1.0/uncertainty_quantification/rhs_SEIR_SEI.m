% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
% Initial value problem: rhs_SEIR_SEI.m
%
% This function defines the system of ODEs for the 
% SEIR-SEI epidemiological model for ZikaV dynamics
%
% The dynamic state coordinates are:
%
%   SH = susceptibles humans          (number of individuals)
%   EH = exposed humans               (number of individuals)
%   IH = infectious humans            (number of individuals)
%   RH = recovered humans             (number of individuals)
%
%   SV = proportion of susceptible vectors    (dimensionless)
%   EV = proportion of exposed vectors        (dimensionless)
%   IV = proportion of infectious vectors     (dimensionless)
%
%   C = cumulative number of infectious humans (number of individuals)
%
% The epidemiological model parameters are:
%
%   N      = human population size (number of individuals)
%
%   betaH  = vector-to-human transmission rate   (days^-1)
%   alphaH = human latent rate                   (days^-1)
%   gamma  = human recovery rate                 (days^-1)
%
%   betaV  = human-to-vector transmission rate   (days^-1)
%   alphaV = vector latent rate                  (days^-1)
%   delta  = vector birth/mortality rate         (days^-1)
%
% Inputs:
%   t: time                    - double
%   y: state vector            - double array (8x1)
%   param: parameters vector   - double array (7x1)
%
% Output:
%   dydt: state rate of change - double array (8x1)
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
% number of lines: 38
% last update: Jun 4, 2021
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function dydt = rhs_SEIR_SEI(t,y,param)

if length(param) < 7 
   error('Warning: To few model parameters')
elseif length(param) > 7
   error('Warning: To many model parameters')
end

if mod(param(1),1) ~= 0
   error('Warning: Use a integer human population size')
end

if sum(param<0)~=0
   error('Warning: Use positive model parameters values')
end

if length(y) < 8 
   error('Warning: To few model coordinates')
elseif length(y) > 8
   error('Warning: To many model coordinates')
end

%if sum(y<0)~=0
%   error('Warning: Use positive model coordinates values')
%end

% model parameters: 
% param = [N betaH alphaH gamma betaV alphaV delta]

N      = param(1);    % human population size (number of individuals)
betaH  = param(2);    % vector-to-human transmission rate (days^-1)
alphaH = param(3);    % human latent rate (days^-1)
gamma  = param(4);    % human recovery rate (days^-1)
betaV  = param(5);    % human-to-vector transmission rate (ind^-1days^-1)
alphaV = param(6);    % incubation vector frequency (days^-1)
delta  = param(7);    % vector lifespan frequency (days^-1)


% SEIR-SEI dynamic model:
%
%      y = [SH EH IH RH SV EV IV C]                         is the state vector
%   dydt = [dSHdt dEHdt dIHdt dRHdt dSVdt dEVdt dIVdt dCdt] is the evolution law
% 
% dSHdt - rate of susceptible humans                (number of individuals/days)
% dEHdt - rate of exposed humans                    (number of individuals/days)
% dIHdt - rate of infectious humans                 (number of individuals/days)
% dRHdt - rate of recovered humans                  (number of individuals/days)
%
% dSHdt - rate of susceptible vectors               (proportion of vectors/days)
% dEHdt - rate of exposed vectors                   (proportion of vectors/days)
% dIHdt - rate of infectious vectors                (proportion of vectors/days)
%
% dCdt - rate of cumulative infectious humans       (number of individuals/days)

[SH EH IH RH SV EV IV C] = deal(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8));

dSHdt = - betaH*SH.*IV;
dEHdt = betaH*SH.*IV - alphaH*EH;
dIHdt = alphaH*EH - gamma*IH;
dRHdt = gamma*IH;
dSVdt = delta - betaV.*SV.*(IH/N) - delta*SV;
dEVdt = betaV*SV.*(IH/N) - alphaV*EV - delta*EV;
dIVdt = alphaV*EV - delta*IV;
dCdt  = alphaH*EH;
dydt = [dSHdt; dEHdt; dIHdt; dRHdt; dSVdt; dEVdt; dIVdt; dCdt];
                                                                        
end
% -----------------------------------------------------------------

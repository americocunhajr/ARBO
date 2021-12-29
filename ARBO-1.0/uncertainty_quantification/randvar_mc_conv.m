
% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
%  randvar_mc_conv.m
%
%  This function computes the metric used to evaluate the 
%  convergence of Monte Carlo method.
%
%              _                                      _ 1/2
%             |                     -- t1              |
%             |  1                 |                   |
%  MC_conv := | --- sum_{j=1}^{Ns} | ||Q(t,s_j)||^2 dt |
%             | Ns                 |                   |
%             |_                 -- t = t0            _|
%  
%  Input:
%  normL2_Q - (1 x Ns) L2 norm of the 2-norm of Q(t)
%  Ns       - number of MC realizations
%
%  Output:
%  MC_conv - (1 x Ns) MC convergence metric
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jan 17, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function MC_conv = randvar_mc_conv(normL2_Q,Ns)
    
    % check number of arguments
    if nargin < 1
        error(' Too few inputs.')
    elseif nargin > 2
        error(' Too many inputs.')
    end
    
    
    % compute the MC convergence metric
    MC_conv = sqrt(cumsum(normL2_Q)./(1:Ns));
    
return
% -----------------------------------------------------------------

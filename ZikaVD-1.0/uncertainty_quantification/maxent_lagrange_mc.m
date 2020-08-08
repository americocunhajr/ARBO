
% -----------------------------------------------------------------
%  maxent_lagrange_mc.m
%
%  This function computes the Lagrange multipliers for a
%  maximum entropy (MaxEnt) distribution given the first
%  N statistical moments constraints.
%
%  Reference:
%  A Matlab Program to Calculate the Maximum Entropy Distributions
%  Author: Ali Mohammad-Djafari
%  In: Smith C.R., Erickson G.J., Neudorfer P.O. (eds)
%  Maximum Entropy and Bayesian Methods, pp 221-233
%  Springer, 1992
%
%  input:
%  xmin    - random variable support left extreme
%  xmin    - random variable support right extreme
%  Nx      - number of points for support discretization
%  mu      - (N x 1) statistical moments
%  lambda0 - (N x 1) Lagrange multipliers initial guess (optional)
%  tol     - numerical tolerance (optional)
%  maxiter - maximum number of iterations (optional)
%
%  output:
%  lambda  - (N  x 1) Lagrange multipliers vector
%  pdf     - (Nx x 1) random variable PDF
%  entropy - (N  x 1) entropy value of the PDF
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jan 21, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [lambda,pdf,supp,entropy] = ...
             maxent_lagrange_mc(xmin,xmax,Nx,mu,lambda0,tol,maxiter)

    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 7
        error('Too many inputs.')
    end
    
    % check arguments
    if xmin >= xmax
        error('xmin must be less than xmax.')
    end
    
    if Nx < 2
        error('Nx must be an integer greater than one.')
    end
    
    % number of constraints
    N = length(mu);
    
    if nargin == 4
        % prealocate memory for lambda0
        lambda0 = zeros(N,1);
        % initial guess for lambda
        lambda0(1) = log(xmax-xmin);
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 5
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 6
        % maximum of iteration
        maxiter = 20;
    end
    
    % check arguments
    if length(mu) ~= length(lambda0)
        error('mu and lambda0 vectors must be same length')
    end
    
    % prealocate memory for entropyX
    entropy = zeros(maxiter,1);
    
    % prealocate memory for phi
	phi = zeros(Nx,2*N-1);
    
    % prealocate memory for G
    G = zeros(2*N-1,1);
    
    % prealocate memory for gnk
    gnk = zeros(N,N);
    
    % discretization of random variable support
    supp = linspace(xmin,xmax,Nx)';
    dx   = supp(2) - supp(1);
    
    % define constraint functions
    % (phi_0(x) = 1 and phi_n(x) = x.^n, n = 1...N)
    phi(:,1) = ones(Nx,1);
    for n = 2:2*N-1
        phi(:,n) = supp.*phi(:,n-1);
    end

    % initialize iteration counter
    iter = 0;
    
    % initial guess for lambda
    lambda = lambda0;

    % Newton method iteration
    while iter < maxiter
        
        % update iteration counter
        iter = iter + 1;
        
        % compute MaxEnt PDF
        pdf = exp(-(phi(:,1:N)*lambda));
        
        % compute nonlinear equations G_n(lambda) = mu_n
        for n = 1:2*N-1
            G(n) = dx*sum(phi(:,n).*pdf);
        end
        
        % compute entropy value
        entropy(iter) = lambda'*G(1:N);
        
        % compute Hankel matrix
        for i = 1:N
            gnk(:,i) = -G(i:N+i-1);
        end

        % compute v
        v = mu - G(1:N);
        
        % compute delta
        delta = gnk\v;
        
        % compute lambda
        lambda = lambda + delta;
        
        % check convergence
        if abs(delta./lambda) < tol
            break
        end
        
        % check convergence
        if iter > 2
            if abs((entropy(iter)-entropy(iter-1))/entropy(iter)) < tol
                break
            end
        end
        
    end
	
	% compute the final PDF
    pdf = exp(-(phi(:,1:N)*lambda));
    
end

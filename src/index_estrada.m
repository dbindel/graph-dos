% EE = index_estrada(c)
% EE = index_estrada(Afun, n, nZ, N)
% EE = index_estrada(A, nZ, N)
%
% Compute the Estrada index tr(exp(A)) by the Chebyshev expansion of 
% exponential function and moments of A; the spectrum of A should 
% already lie in [-1,1]. 
%
% Inputs:
%    c: A column vector of N Chebyshev moment estimates of A
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%
% Output:
%    EE: Estrada index of the form tr(exp(A))
%
function EE = index_estrada(varargin)
    % Use the moments if they are given; otherwise compute them
    if min(size(varargin{1}))==1 && length(varargin{1})>1
        c = varargin{1};
    else
        [c,~] = moments_cheb_dos(varargin{:});
    end
    
    % Get the coefficients of Chebyshev expansion of exponential
    N = length(c);
    w = moments_exponential(N);
    EE=w'*c;
    
end
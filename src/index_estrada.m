% EE = index_estrada(c, ab )
% EE = index_estrada(Afuns, n, nZ, N, ab)
% EE = index_estrada(As, nZ, N, ab)
%
% Compute the Estrada index tr(exp(A)) by the Chebyshev expansion of 
% exponential function and moments of A; the spectrum of As should 
% already lie in [-1,1], and the rescaling parameters should be given
% as the last argument.
%
% Inputs:
%    c: A column vector of N Chebyshev moment estimates of A
%    As: Rescaled matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    ab: Rescaling parameters. Original matrix should be ab(1)*A+ab(2)
%
% Output:
%    EE: Estrada index of the form tr(exp(A))
%
function EE = index_estrada(varargin)
    % If rescaling parameters are given
    if length(varargin{end})==2
        ab = varargin{end};
        varargin = varargin(1:end-1);
    else
        ab = [1,0];
    end
    
    % Use the moments if they are given; otherwise compute them
    if min(size(varargin{1}))==1 && length(varargin{1})>1
        c = varargin{1};
    else
        [c,~] = moments_cheb_dos(varargin{:});
    end
    
    % Get the coefficients of Chebyshev expansion of exponential
    N = length(c);
    w = moments_exponential(N,ab(1));
    EE=exp(ab(2))*w'*c;
    
end
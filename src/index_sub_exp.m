% ESC = index_sub_exp(c)
% ESC = index_sub_exp(beta, c)
% ESC = index_sub_exp(beta, Afun, n, nZ, N)
% ESC = index_sub_exp(beta, A, nZ, N)
%
% Compute the exponential subgraph centrality diag(exp(beta*A)) by the
% Chebyshev expansion of exponential function and moments of A; the 
% spectrum of A should already lie in [-1,1]. 
%
% Inputs:
%    beta: Exponential coefficient
%    c: An N-by-n matrix of ldos moment estimates of A
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%
% Output:
%    ESC: Exponential subgraph centraility of the form diag(exp(beta*A))
%
function ESC = index_sub_exp(varargin)
    % Check if an exponential coefficient is given (default:1)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        beta = varargin{1};
        varargin = varargin(2:end);
    else
        beta = 1;
    end
    
    % Check if ldos moments are given; otherwise calculate an estimate
    if size(varargin{1},1)~=size(varargin{1},2)
        c = varargin{1};
    else
        [c,~] = moments_cheb_ldos(varargin{:});
    end
    
    % Get the Chebyshev expansion of exp(beta*x)
    N = size(c,1);
    w = moments_exponential(N,beta);
    ESC = c'*w;
    
end
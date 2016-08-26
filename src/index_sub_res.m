% ESC = index_sub_res(c)
% ESC = index_sub_res(alpha, c)
% ESC = index_sub_res(alpha, Afun, n, nZ, N)
% ESC = index_sub_res(alpha, A, nZ, N)
%
% Compute the resolvent subgraph centrality diag((I-alpha*A)^-1) by 
% the Chebyshev expansion of exponential function and moments of A; 
% the spectrum of A should already lie in [-1,1]. 
%
% Inputs:
%    alpha: Resolvent coefficient
%    c: An N-by-n matrix of ldos moment estimates of A
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%
% Output:
%    RSC: Resolvent subgraph centraility of the form diag((I-alpha*A)^-1)
%
function RSC = index_sub_res(varargin)
    % Check if an exponential coefficient is given (default:0.9)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        alpha = varargin{1};
        varargin = varargin(2:end);
    else
        alpha = 0.9;
    end
    
    % Check if ldos moments are given; otherwise calculate an estimate
    if size(varargin{1},1)~=size(varargin{1},2)
        c = varargin{1};
    else
        [c,~] = moments_cheb_ldos(varargin{:});
    end
    
    % Get the Chebyshev expansion of 1/(1-alpha*x)
    N = size(c,1);
    w = moments_resolvent(N,alpha);
    RSC = c'*w;
    
end
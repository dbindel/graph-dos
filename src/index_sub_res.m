% ESC = index_sub_res(c, ab)
% ESC = index_sub_res(alpha, c, ab)
% ESC = index_sub_res(alpha, Afuns, n, nZ, N, ab)
% ESC = index_sub_res(alpha, As, nZ, N, ab)
%
% Compute the resolvent subgraph centrality diag((I-alpha*A)^-1) by 
% the Chebyshev expansion of exponential function and moments of A; 
% the spectrum of As should already lie in [-1,1], and the rescaling 
% parameters should be given as the last argument. 
%
% Inputs:
%    alpha: Resolvent coefficient
%    c: An N-by-n matrix of ldos moment estimates of A
%    As: Rescaled matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    ab: Rescaling parameters. Original matrix should be ab(1)*A+ab(2)
%
% Output:
%    RSC: Resolvent subgraph centraility of the form diag((I-alpha*A)^-1)
%
function RSC = index_sub_res(varargin)
    % If rescaling parameters are given
    if length(varargin{end})==2
        ab = varargin{end};
        varargin = varargin(1:end-1);
    else
        ab = [1,0];
    end
    
    % Check if an exponential coefficient is given (default:0.9)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        alpha = varargin{1};
        varargin = varargin(2:end);
    else
        alpha = 0.9;
    end
    
    % Throw a warning if choice of alpha does not guarantee convergence
    if alpha>=1/sum(ab)
        warning(['Chebyshev expansion will not converge because'...
                ' function has pole in [-1,1]. alpha must be less'...
                ' than 1/sum(ab)']);
    end
    
    % Check if ldos moments are given; otherwise calculate an estimate
    if size(varargin{1},1)~=size(varargin{1},2)
        c = varargin{1};
    else
        [c,~] = moments_cheb_ldos(varargin{:});
    end
    
    % Get the Chebyshev expansion of 1/(1-alpha*x)
    N = size(c,1);
    w = moments_resolvent(N,alpha,ab);
    RSC = c'*w;
    
end
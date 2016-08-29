% ESC = index_sub_exp(c, ab)
% ESC = index_sub_exp(beta, c, ab)
% ESC = index_sub_exp(beta, Afuns, n, nZ, N, ab)
% ESC = index_sub_exp(beta, As, nZ, N, ab)
%
% Compute the exponential subgraph centrality diag(exp(beta*A)) by the
% Chebyshev expansion of exponential function and moments of A; the 
% spectrum of As should already lie in [-1,1], and the rescaling 
% parameters should be given as the last argument.
%
% Inputs:
%    beta: Exponential coefficient
%    c: An N-by-n matrix of ldos moment estimates of A
%    As: Rescaled matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    ab: Rescaling parameters. Original matrix should be ab(1)*A+ab(2)
%
% Output:
%    ESC: Exponential subgraph centraility of the form diag(exp(beta*A))
%
function ESC = index_sub_exp(varargin)
    % If rescaling parameters are given
    if length(varargin{end})==2
        ab = varargin{end};
        varargin = varargin(1:end-1);
    else
        ab = [1,0];
    end
    
    % Check if an exponential coefficient is given (default:1)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        beta = varargin{1};
        varargin = varargin(2:end);
    else
        beta = 1;
    end
    
    ab = beta*ab;
    
    % Check if ldos moments are given; otherwise calculate an estimate
    if size(varargin{1},1)~=size(varargin{1},2)
        c = varargin{1};
    else
        [c,~] = moments_cheb_ldos(varargin{:});
    end
    
    % Get the Chebyshev expansi0n of exp(beta*x)
    N = size(c,1);
    w = moments_exponential(N,ab(1));
    ESC = exp(ab(2))*(c'*w);
    
end
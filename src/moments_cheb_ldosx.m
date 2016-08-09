% [c, cs] = moments_cheb_ldosx(Afun, n, nodes, N, kind)
% [c, cs] = moments_cheb_ldosx(A, nodes, N, kind)
%
% Compute a column vector (or vectors) of Chebyshev moments of
% the form c(k,j) = [T_k(A)]_jj for k in nodes.  This routine
% does no scaling; the spectrum of A should already lie in [-1,1].
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nodes: Indides of entries to compute
%    N: Number of moments to compute
%    kind: 1 or 2 for first or second kind Chebyshev functions
%          (default is 1)
%
% Output:
%    c: an N-by-n matrix of moments
%
function [c] = moments_cheb_ldosx(A, n, nodes, N, kind)

if isa(A, 'function_handle')
  % Fill in defaults
  Afun = A;
  if nargin < 5, kind = 1;    end
  if nargin < 4, N  = 10;     end
  if nargin < 3, nodes = [1]; end
  if N < 2,      N = 2;       end
else
  % Fill in defaults and shift arguments down
  Afun = @(X) A*X;
  if nargin < 4, kind = 1;    else kind = N;    end
  if nargin < 3, N = 10;      else N = nodes;   end
  if nargin < 2, nodes = [1]; else nodes = n;   end
  n = size(A,1);
end

% Set up random probe vectors (we allow them to be passed in, too)
V = zeros(n, length(nodes));
V(nodes,:) = eye(length(nodes));

% Estimate moments for each probe vector
c = moments_cheb(Afun, V, N, kind);

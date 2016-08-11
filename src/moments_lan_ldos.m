% [c] = moments_lan_ldos(Afun, n, nodes, N, kmax, btol)
% [c] = moments_lan_ldos(A, nodes, N, kmax, btol)
%
% Compute a column vector (or vectors) of Chebyshev moments of
% the form c(k) = tr(T_k(A)) for k = 0 to N-1.  This routine
% does no scaling; the spectrum of A should already lie in [-1,1].
% The traces are computed via a stochastic estimator with nZ probes.
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    N: Number of moments to compute
%    kmax: maximum Lanczos steps (default is 100)
%    btol: Tolerance on Lanczos residual (default is 1e-6)
%
% Output:
%    c:  Chebyshev moments
%
function [c] = moments_lan_ldos(A, n, nodes, N, kmax, btol)

  if isa(A, 'function_handle')
    % Fill in defaults
    Afun = A;
    if nargin < 6, btol  = 1e-6; end
    if nargin < 5, kmax  = 100;  end
    if nargin < 4, N     = 10;   end
    if nargin < 3, nodes = [1];  end
  else
    % Fill in defaults and shift arguments down
    Afun = @(X) A*X;
    if nargin < 5, btol  = 1e-6; else btol  = kmax;  end
    if nargin < 4, kmax  = 100;  else kmax  = N;     end
    if nargin < 3, N     = 10;   else N     = nodes; end
    if nargin < 2, nodes = [1];  else nodes = n;     end
    n = size(A,1);
  end

  if length(nodes) < n
    V = zeros(n,length(nodes));
    V(nodes,:) = eye(length(nodes));
  else
    V = nodes;
  end

  c = zeros(N, size(V,2));
  for k = 1:size(nodes,2)
    [theta,wts] = moments_lanczos(Afun, n, V(:,k), kmax, btol);
    c(:,k) = moments_quad2cheb(theta, wts, N);
  end

end

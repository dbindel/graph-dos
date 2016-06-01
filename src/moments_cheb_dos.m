% [c, cs] = moments_cheb_dos(Afun, n, nZ, N, kind)
% [c, cs] = moments_cheb_dos(A, nZ, N, kind)
%
% Compute a column vector (or vectors) of Chebyshev moments of
% the form c(k) = tr(T_k(A)) for k = 0 to N-1.  This routine
% does no scaling; the spectrum of A should already lie in [-1,1].
% The traces are computed via a stochastic estimator with nZ probes.
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    kind: 1 or 2 for first or second kind Chebyshev functions
%          (default is 1)
%
% Output:
%    c: an column vector of N moment estimates
%    cs: standard deviation of the moment estimator (std/sqrt(nZ))
%
function [c,cs] = moments_cheb_dos(A, n, nZ, N, kind)

  if isa(A, 'function_handle')
    % Fill in defaults
    Afun = A;
    if nargin < 5, kind = 1;  end
    if nargin < 4, N  = 10;   end
    if nargin < 3, nZ = 100;  end
    if N < 2,      N = 2;     end
  else
    % Fill in defaults and shift arguments down
    Afun = @(X) A*X;
    if nargin < 4, kind = 1; else kind = N; end
    if nargin < 3, N = 10;   else N = nZ;   end
    if nargin < 2, nZ = 100; else nZ = n;   end
    n = size(A,1);
  end

  % Set up random probe vectors (we allow them to be passed in, too)
  if length(nZ) > 1
    Z = nZ;
    nZ = size(Z,2);
  else
    Z = randn(n,nZ);
  end

  % Estimate moments for each probe vector
  cZ = moments_cheb(A, Z, N, kind);
  c  = mean(cZ,2);
  if nargout > 1
    cs = std(cZ,0,2)/sqrt(nZ);
  end
  
end

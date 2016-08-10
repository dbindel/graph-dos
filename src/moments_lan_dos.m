% [theta, wts] = moments_lan_dos(Afun, n, nZ, N, kmax, btol)
% [theta, wts] = moments_lan_dos(A, nZ, N, kmax, btol)
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
%    kmax: maximum Lanczos steps (default is 100)
%    btol: Tolerance on Lanczos residual (default is 1e-6)
%
% Output:
%    c:  Chebyshev moments
%    cs: Standard deviations
%
function [c, cs] = moments_lan_dos(A, n, nZ, N, kmax, btol)

  if isa(A, 'function_handle')
    % Fill in defaults
    Afun = A;
    if nargin < 6, btol = 1e-6;  end
    if nargin < 5, kmax = 100;   end
    if nargin < 4, N    = 10;    end
    if nargin < 3, nZ   = 100;   end
  else
    % Fill in defaults and shift arguments down
    Afun = @(X) A*X;
    if nargin < 5, btol = 1e-6; else btol = kmax; end
    if nargin < 4, kmax = 100;  else kmax = N;    end
    if nargin < 3, N    = 10;   else N    = nZ;   end
    if nargin < 2, nZ   = 100;  else nZ   = n;    end
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
  cZ = zeros(N,nZ);
  for k = 1:nZ
    [theta,wts] = moments_lanczos(Afun, n, Z(:,k), kmax, btol);
    cZ(:,k) = moments_quad2cheb(theta, wts, N);
  end

  % Get summary statistics
  c = mean(cZ,2);
  if nargout > 1
    cs = std(cZ,0,2)/sqrt(nZ);
  end

end

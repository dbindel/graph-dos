% [c, cs] = moments_cheb(A, V, N, kind)
%
% Compute a column vector of Chebyshev moments of
% the form c(k) = v'*T_k(A)*v for k = 0 to N-1.  This routine
% does no scaling; the spectrum of A should already lie in [-1,1].
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    V: Starting vectors
%    N: Number of moments to compute
%    kind: 1 or 2 for first or second kind Chebyshev functions
%          (default is 1)
%
% Output:
%    c: a length N vector of moments
%
function [c] = moments_cheb(A, V, N, kind)

  % Set parameter defaults; we require at least two moments
  if nargin < 4, kind = 1;  end
  if nargin < 2, N = 10;    end
  if N < 2,      N = 2;     end

  % We only need matvecs
  if isa(A, 'function_handle')
    Afun = A;
  else
    Afun = @(X) A*X;
  end

  [n,p] = size(V);
  c = zeros(N,p);

  % Run three-term recurrence to compute moments
  TVp = V;
  TVk = kind*Afun(V);
  c(1,:) = sum(V.*TVp);
  c(2,:) = sum(V.*TVk);
  for k = 3:N
    TV = 2*Afun(TVk) - TVp;
    TVp = TVk;
    TVk = TV;
    c(k,:) = sum(V.*TVk);
  end

end

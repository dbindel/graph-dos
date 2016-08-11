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
function [c] = moments_cheb_ldosx(varargin)

  defaults = {'Afun', NaN, 'n', NaN, ...
              'nodes', NaN, 'N', 10, 'kind', 1};
  [Afun, n, nodes, N, kind] = mfuncify(defaults, varargin{:});
  if N < 2, N = 2; end

  % Set up random probe vectors (we allow them to be passed in, too)
  V = zeros(n, length(nodes));
  V(nodes,:) = eye(length(nodes));

  % Estimate moments for each probe vector
  c = moments_cheb(Afun, V, N, kind);

end

% [Lfun] = mfunc_slaplacian(Wfun, n)
%
% Convert a weighted adjacency matrix function into a signless Laplacian.
%
% Input:
%   Wfun: multiplies by a weighted adjacency matrix
%   n: dimension of the space
%
% Output:
%   Lfun: multiplies by a graph Laplacian
%
function [Lfun] = mfunc_slaplacian(W, n)

  if isa(W, 'function_handle')
    Wfun = W;
  else
    n = size(W,1);
    Wfun = @(X) W*X;
  end

  e = ones(n,1);
  D = spdiags(Wfun(e), 0, n, n);
  Lfun = @(X) (D*X + Wfun(X));

end

% [Bfun] = mfunc_modularity(Wfun, n)
%
% Convert a weighted adjacency matrix function into a modularity matrix.
%
% Input:
%   Wfun: multiplies by a weighted adjacency matrix
%   n: dimension of the space
%
% Output:
%   Bfun: multiplies by a graph modularity matrix
%
function [Bfun] = mfunc_modularity(W, n)

  if isa(W, 'function_handle')
    Wfun = W;
  else
    n = size(W,1);
    Wfun = @(X) W*X;
  end

  e = ones(n,1);
  d = Wfun(e);
  s = sum(d);
  Bfun = @(X) (Wfun(X) - (d/s)*(d'*X));

end

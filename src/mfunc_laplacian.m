% [Lfun] = mfunc_laplacian(Wfun, n)
%
% Convert a weighted adjacency matrix function into a Laplacian.
%
% Input:
%   Wfun: multiplies by a weighted adjacency matrix
%   n: dimension of the space
%
% Output:
%   Lfun: multiplies by a graph Laplacian
%
function [Lfun] = mfunc_laplacian(Wfun, n)

  if ~isa(Wfun, 'function_handle')
    n = size(Wfun,1);
    Wfun = @(X) W*X;
  end

  e = ones(n,1);
  D = spdiags(Wfun(e), 0, n, n);
  Lfun = @(X) (D*X - Wfun(x));

end

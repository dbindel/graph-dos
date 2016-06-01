% [Nfun] = mfunc_normalize(Wfun, n, mode)
%
% Normalize a weighted adjacency matrix.  Note that for anything
% but row normalization, we assume the matrix is symmetric.
%
% Input:
%   W: function to apply weighted adjacency matrix
%   n: dimension of the space
%   mode: string indicating the style of normalization; this may be
%     's': Symmetric scaling by the degrees (default)
%     'r': Normalize to row-stochastic
%     'c': Normalize to col-stochastic
%
% Output:
%   Nfun: applies normalized adjacency matrix or stochastic matrix
%
function [Nfun] = mfunc_adjacency(Wfun, n, mode)

  if ~isa(Wfun, 'function_handle')
    n = size(Wfun,1);
    Wfun = @(X) W*X;
  end
  if nargin < 3, mode = 's'; end

  e = ones(n,1);
  d = Wfun(e);

  if strcmp(mode, 's')
    D = spdiags(1./sqrt(d), 0, n, n);
    Nfun = @(X) D*Wfun(D*X);
  elseif strcmp(mode, 'r')
    D = spdiags(1./d, 0, n, n);
    Nfun = @(X) D*Wfun(X);
  elseif strcmp(mode, 'c')
    D = spdiags(1./d, 0, n, n);
    Nfun = @(X) Wfun(D*X);
  else
    error('Unknown mode!');
  end

end

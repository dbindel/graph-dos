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
function [Nfun] = mfunc_normalize(W, n, mode)

  if isa(W, 'function_handle')
    if nargin < 3, mode = 's'; end
    Wfun = W;
  else
    if nargin < 2
      mode = 's';
    elseif nargin < 3
      if ischar(n), mode = n; else mode = 's'; end
    end
    n = size(W,1);
    Wfun = @(X) W*X;
  end

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

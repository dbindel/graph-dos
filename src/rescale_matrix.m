% [Hs, ab] = rescale_matrix(H, n, range)
% [Hs, ab] = rescale_matrix(H, n, fudge)
%
% Rescale the symmetric matrix H so the eigenvalue range maps to
% between -1 and 1.  If a range is not given, the function uses Lanczos
% to estimate the extremal eigenvalues, expanding by a relative fudge
% factor.
%
% Input:
%    H: The original matrix (or function)
%    n: The dimension of the space (iff H is a function)
%    range: A two-element vector representing an interval where eigs live
%    fudge: A scalar "fudge factor" to guard against range underestimates
%           (default is 0.01)
%
% Output:
%    Hs: A function representing the scaled matrix
%    ab: Transformation parameters: Hs = (H-b)/a
%
function [Hs, ab] = rescale_matrix(H, n, range)

  % Deal with matrix vs function
  if isa(H, 'function_handle')
    Hfun = H;
    if nargin < 2, error('Missing size argument'); end
    if nargin < 3, range = 0.01; end
  else
    Hfun = @(x) H*x;
    n = size(H,1);
    if nargin < 2, range = 0.01; else range = n; end
  end

  % Run Lanczos to estimate range (if needed)
  if length(range) == 1
    fudge = range;
    opts = [];
    opts.isreal = 1;
    opts.issym = 1;
    range = sort(eigs(Hfun, n, 2, 'be', opts));
  else
    fudge = 0;
  end

  % Parameters for linear mapping
  ab = [(range(2)-range(1))/(2-fudge); (range(2)+range(1))/2];

end

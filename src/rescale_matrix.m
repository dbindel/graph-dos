% [H, ab] = rescale_matrix(H, range)
% [H, ab] = rescale_matrix(H, fudge)
%
% Rescale the symmetric matrix H so the eigenvalue range maps to
% between -1 and 1.  If a range is not given, the function uses Lanczos
% to estimate the extremal eigenvalues, expanding by a relative fudge
% factor.
%
% Input:
%    H: The original matrix
%    range: A two-element vector representing an interval where eigs live
%    fudge: A scalar "fudge factor" to guard against range underestimates
%           (default is 0.01)
%
% Output:
%    H: The scaled matrix
%    ab: Transformation parameters: Hs = (H-b)/a
%
function [H, ab] = rescale_matrix(H, range)

  % Get the fudge factor and range
  if nargin < 2
    fudge = 0.01;
    range = [];
  elseif length(range) == 1
    fudge = range;
    range = [];
  else
    fudge = 0;
  end

  % Compute range if not given
  if isempty(range)
    opts = [];
    opts.isreal = 1;
    opts.issym = 1;
    range = sort(eigs(H, 2, 'be', opts));
  end

  % Form dense or sparse identity, as appropriate
  if issparse(H)
    I = speye(length(H));
  else
    I = eye(length(H));
  end

  % Parameters for linear mapping
  ab = [(range(2)-range(1))/(2-fudge); (range(2)+range(1))/2];
  H = (H-ab(2)*I)/ab(1);
end

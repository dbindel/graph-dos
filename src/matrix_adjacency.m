% [A] = matrix_adjacency(W, keep_diag)
%
% Convert a weighted adjacency matrix into an ordinary adjacency.
%
% Input:
%   W: weighted adjacency matrix
%   keep_diag: flag if we want to keep the diagonal on output (default: 1)
%
% Output:
%   A: a 0-1 adjacency matrix (in sparse form)
%
function [A] = matrix_adjacency(W, keep_diag)

  if nargin < 2, keep_diag = 1; end

  % Zero out the diagonal if it is not wanted
  if ~keep_diag
    z = zeros(size(W,1), 1);
    W = spdiags(z, 0, W);
  end

  % Reconstruct with only off-diagonals
  [i,j,wij] = find(W);
  A = sparse(i, j, ones(size(wij)), size(W,1), size(W,2));

end

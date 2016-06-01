% [L] = matrix_laplacian(W, mode)
%
% Convert a weighted adjacency matrix into a Laplacian.
%
% Input:
%   W: weighted adjacency matrix
%   mode: 'r'ow or 'c'ol sum to zero (default is rows)
%
% Output:
%   L: a graph Laplacian
%
function [L] = matrix_laplacian(W, mode)

  if nargin < 2, mode = 'r'; end

  % Strip diagonal
  z = zeros(size(W,1),1);
  W = spdiags(z, 0, W);

  % Compute row or column sums
  if strcmp(mode, 'r')
    d = full(sum(W,2));
  else
    d = full(sum(W,1)).';
  end

  % Replace main diagonal, flip off-diagonal signs
  L = spdiags(d, 0, -W);

end

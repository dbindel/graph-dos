% [xy] = gplot_spectral(A)
%
% Spectral graph layout and plotting.
%
function [xy] = gplot_spectral(A)
  N = matrix_normalize(A);

  % Deal with different sizes
  if length(A) > 100
    [V,D] = eigs(A);
    xy = V(:,2:3);
  elseif length(A) >= 3
    [V,D] = eig(full(A));
    xy = V(:,2:3);
  elseif length(A) == 2
    xy = [-0.5,0; 0.5,0];
  else
    xy = [0,0];
  end

  % Plot if no output
  if nargout > 1
    gplot(A, xy, '-*');
  end
end

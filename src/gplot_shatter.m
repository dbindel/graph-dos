% [xy] = gplot_shatter(A, min_size, xy)
%
% Spectral graph layout and plotting.  Assumes there may be several
% disconnected components, and plots them in a grid for visual
% comparison.
%
% Args:
%   A: Adjacency matrix (potentially scaled)
%   min_size: Smallest size component to include
%   xy: Coordinates (optional, use spectral coords otherwise)
%
% Output:
%   xy: Rewritten coordinates; nodes for small components map to zeros
%
function [xy] = gplot_shatter(A, min_size, xy)

  if nargin < 2, min_size = 2;             end
  if nargin < 3, xy = zeros(length(A), 2); end

  [nc,sizes,members] = networkComponents(A);
  nc_big = sum(sizes >= min_size);
  nrow = ceil(sqrt(nc_big));

  jj = 1;
  for j = 1:nc
    Ij = members{j};
    col = mod(jj-1,nrow);
    row = floor((jj-1)/nrow);
    jj = jj+1;

    if sizes(j) >= min_size

      % Get coordinates for component
      if nargin < 3
        xyj = gplot_spectral(A(Ij,Ij));
      else
        xyj = xy(Ij,:);
      end

      % Rescale to fit in [-0.4,0.4]^2 + [col,row]
      xyrange = 1e-8 + max(xyj)-min(xyj);
      xymean  = mean(xyj);
      xyj(:,1) = col + 0.5 + 0.4*(xyj(:,1)-xymean(1))/xyrange(1);
      xyj(:,2) = row + 0.5 + 0.4*(xyj(:,2)-xymean(2))/xyrange(2);

      % Write back coordinates
      xy(Ij,:) = xyj;

    else
      xy(Ij,:) = 0;
    end
  end

  % Plot if no output
  if nargout < 1
    fprintf(' %d', sizes(sizes >= min_size));
    fprintf('\nToo small: %d\n', nc-nc_big);
    gplot(A, xy, '-*');
  end
end

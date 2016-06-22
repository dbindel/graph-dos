% [c] = moments_delta(x,N)
%
% Compute Chebyshev moments 0 through N-1 of a delta function
% centered at x.
%
function [c] = moments_delta(x,N)
  c = cos( (0:N-1).*acos(x) )';
end

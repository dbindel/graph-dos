% yy = plot_chebint(c,xx,ab)
%
% Given a (possibly filtered) set of first-kind Chebyshev moments,
% compute the integral of the density
%
%   int_{0}^s (2/pi)*sqrt(1-x^2)) * ( c(0)/2 + sum_{n=1}^{N-1} c_n T_n(x) )
%
% If no output argument is assigned, make a plot.

function yy = plot_chebint(c,xx0,ab)

  % Map points to [-1,1]
  if nargin < 3
    ab = [1, 0];
    xx = xx0;
  else
    xx = (xx0-ab(2))/ab(1);
  end

  N = length(c);
  txx = acos(xx);
  yy = c(1)*(txx-pi)/2;
  for np = 2:N
    n = np-1;
    yy = yy + c(np) * sin(n*txx)/n;
  end
  yy = -2/pi * yy;

  if nargout < 1
    plot(xx0, yy);
  end

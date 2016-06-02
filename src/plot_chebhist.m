% yy = plot_chebhist(c,xx)
%
% Given a (possibly filtered) set of first-kind Chebyshev moments,
% compute the integral of the density
%
%   int_{0}^s (2/pi)*sqrt(1-x^2)) * ( c(0)/2 + sum_{n=1}^{N-1} c_n T_n(x) )
%
% If no output argument is assigned, make a plot.

function yy = plot_chebhist(c,xx,varargin)

  yy = plot_chebint(c,xx,varargin{:});
  yy = yy(2:end)-yy(1:end-1);
  xm = (xx(2:end)+xx(1:end-1))/2;
  if nargout < 1
    bar(xm, yy);
  end

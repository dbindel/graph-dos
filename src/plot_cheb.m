% yy = plot_cheb(c,xx,ab)
%
% Given a set of fist-kind Chebyshev moments, compute the associated
% density.   If no output argument is assigned, make a plot.
%
% Inputs:
%   c:  Chebyshev moments (on [-1,1])
%   xx: evaluation points (defaults to mesh of 1001 pts)
%   ab: mapping parameters (defaults to identity)
%
% Output:
%   yy: Density evaluated at xx mesh

function yy = plot_cheb(c,xx0,ab,kind)

  % In practice, allow kind = 2 -- though things like the Jackson
  % filter are set up only for first-kind moments, so we will not
  % advertise this capability in the public comments.
  %
  if nargin < 4, kind = 1;                            end
  if nargin < 2, xx0 = linspace(-1-1e-8,1+1e-8,1001); end

  % Map points to [-1,1]
  if nargin < 3
    ab = [1, 0];
    xx = xx0;
  else
    xx = (xx0-ab(2))/ab(1);
  end

  % Run the recurrence
  N = length(c);
  P0 = xx*0+1;
  P1 = kind*xx;
  yy = c(1)/(3-kind) + c(2)*xx;
  for np = 3:N
    Pn = 2*(xx.*P1) - P0;
    yy = yy + c(np)*Pn;
    P0 = P1;
    P1 = Pn;
  end

  % Do the normalization
  if kind == 1
    yy = (2/pi/ab(1))*( yy./(1e-12+sqrt(1-xx.^2)) );
  else
    yy = (2/pi/ab(1))*( yy.*(sqrt(1-xx.^2)) );
  end

  % Plot if appropriate
  if nargout < 1
    plot(xx0, yy);
  end

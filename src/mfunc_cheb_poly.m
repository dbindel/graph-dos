% [pWfun] = mfunc_cheb_poly(coeff, W)
%
% Given a matrix W (or a function that multiplies by a matrix) and
% a vector of coefficients, return a function that evaluates
%    p(W)*v = sum_{j=1}^d coeff(j) * T_{j-1}(W) * v
%
function [pWfun] = mfunc_cheb_poly(coeff, W)

  if isa(W, 'function_handle')
    Wfun = W;
  else
    Wfun = @(X) W*X;
  end

  if length(coeff) == 1
    pWfun = @(x) coeff * Wfun(x)
  else
    pWfun = @(x) mfunc_cheb_poly_apply(coeff, Wfun, x);
  end

end

function [y] = mfunc_cheb_poly_apply(coeff, Wfun, x)
  Pm = x;
  P0 = Wfun(x);
  y = coeff(1)*Pm;
  for j = 2:length(coeff)-1
    y = y + coeff(j) * P0;
    Pn = 2*Wfun(P0) - Pm;
    Pm = P0;
    P0 = Pn;
  end
  y = y + coeff(end) * P0;
end

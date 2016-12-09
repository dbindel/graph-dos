% Given a matrix {A}, an eigenvalue {lambda}, and
% the support {s} of the corresponding eigenbasis V,
% return V.
function [L,V] = eigenbasis(A, lambda, s, rtol)

  % Form block matrix reordering rows
  I = find(s);
  
  [Lambda, Q] = eig(full(A(I, I)));
  
  M = A(:,I);
  V = zeros(length(A), 1);
  L = zeros(1, 1);
  
  for j=1:length(Lambda)
      l = Lambda(j,j);
      if (abs(lambda - l) < rtol)      
          v = zeros(length(A),1);
          v(I) = Q(:,j);
          resid = norm (M * v(I) - l * v);
          if (resid < rtol)
              L = [L; l];
              V = [V; v];
          end
      end
  end
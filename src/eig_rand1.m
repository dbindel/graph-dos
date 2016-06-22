% [V,D] = eig_rand1(Z, AZ, thresh)
% [V,D] = eig_rand1(Z, Afun, thresh)
%
% If A is an (approximately) low-rank symmetric matrix, estimate eigenpairs
% of A associated from multiplication by a random probe matrix.  Uses a
% one-pass approach -- this could be more accurate with repeated multiplication
% by A.
%
% Args:
%   Z: A random probe matrix
%   AZ: A*Z
%   Afun: A function used to apply A to a block of vectors
%   thresh: Singular value decay threshold used to truncate A*Z
%
% Ref: Uses algorithm 5.6 from Halko, Martinsson, and Tropp.
% (http://arxiv.org/pdf/0909.4061v2.pdf)

function [V,D] = eig_rand1(Z, AZ, thresh)

  % Default to a very tight threshold
  if nargin < 3, thresh = 1e-8; end

  % Apply
  if isa(AZ, 'function_handle'), AZ = AZ(Z); end

  % Get Q an approximate basis for the range space
  [U,S,V] = svd(AZ,0);
  s = diag(S);
  m = sum(s > thresh*s(1));
  Q = U(:,1:m);

  % Should have B*(Q'*Z) = Q'*AZ; choose B to minimize least sq diff
  % The stationary equations look like
  %    symm(B*M*M'-N*M') = 0
  % or
  %    B*M*M' + M*M'*B = N*M' + M*N'
  % where M = Q'*Z and N = Q'*AZ.  This is a Sylvester equation,
  % which we solve via a Bartels-Stewart approach.

  M = Q'*Z;
  N = Q'*AZ;
  C = N*M';
  C = C+C';

  [U,S,V] = svd(M',0);
  s2 = diag(S).^2;
  e = ones(size(s2));
  Bt = (V'*C*V) ./ (s2 * e' + e * s2');
  B = V*Bt*V';

  % Reconstruct the eigenpairs
  [V,D] = eig(B);
  [~,I] = sort(abs(diag(D)), 'descend');
  D = D(I,I);
  V = Q*V(:,I);

end

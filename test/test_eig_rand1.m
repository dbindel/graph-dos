% Form low-rank test matrix
Acore = randn(10,10);
Acore = (Acore+Acore')/2;
[U,R] = qr(randn(1000,10),0);
Afun = @(x) U*(Acore*(U'*x));

% Eigenvalues from random probe
Z = randn(1000,10);
AZ = Afun(Z);
[V,D] = eig_rand1(Z,AZ,1e-6);

% Compute residual
R = Afun(V)-V*D;
fprintf('Residual for recovered eigenvectors: %e\n', norm(R, 'fro'));

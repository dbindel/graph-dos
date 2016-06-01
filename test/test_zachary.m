% Test basic functionality on the Zachary network
clear

% Load adjacency matrix and check basic
A = load_zachary();
m = 78;  % Number of edges
n = 34;  % Number of nodes
assert(nnz(A) == 2*m,  'Wrong number of edges');
assert(size(A,1) == n, 'Wrong number of nodes');

% Compute full eigendecomposition for later sanity checks
[V,D] = eig(full(A));
lambda = diag(D);

% Check symmetric normalization
N = matrix_normalize(A, 's');
lambdaN = eig(full(N));
assert(norm(N-N',1) < 1e-12, 'Symmetric scaling fails to preserve symmetry');
assert(abs(max(lambdaN)-1) < 1e-12, 'Missing eigenvalue near 1 after normal?');
assert(norm(A-matrix_adjacency(N,0),1) == 0, 'Garbled recovery after s scale');

% Check row normalization
P = matrix_normalize(A, 'r');
assert(norm(sum(P,2)-1) < 1e-12, 'Not row normalized');
assert(norm(A-matrix_adjacency(N,0),1) == 0, 'Garbled recovery after r scale');

% Check col normalization
P = matrix_normalize(A, 'c');
assert(norm(sum(P,1)-1) < 1e-12, 'Not col normalized');
assert(norm(A-matrix_adjacency(N,0),1) == 0, 'Garbled recovery after c scale');

% Check adjacency conversion routines
Ab = matrix_normalize(A) + speye(n);
assert(norm(A-matrix_adjacency(Ab,0), 1) == 0, 'Incorrect diagonal discard');
assert(norm(A-matrix_adjacency(Ab,1), 1) == 1, 'Incorrectly discarded diag');

% Check Laplacian conversion routines
L = matrix_laplacian(A);
assert(norm(sum(L,2)) < 1e-12, 'Incorrect row sum');
assert(norm(A-matrix_adjacency(L,0),1) == 0, 'Incorrect Laplacian -> adj');

% Check Laplacian function conversion
Lfun = mfunc_laplacian(A);
Lfun2 = mfunc_laplacian(@(x) A*x, n);
x = randn(n,1);
assert(norm(Lfun(x)-L*x) < 1e-12*norm(L*x), 'Incorrect Lfun');
assert(norm(Lfun2(x)-L*x) < 1e-12*norm(L*x), 'Incorrect Lfun2');

% Check function normalization functionality (row scaling)
mode = {'r', 'c', 's'};
for k = 1:3
  P = matrix_normalize(A, mode{k});
  Pfun = mfunc_normalize(A, mode{k});
  Pfun2 = mfunc_normalize(@(x) A*x, n, mode{k});
  assert(norm(Pfun(x)-P*x) < 1e-12*norm(P*x), 'Incorrect Pfun');
  assert(norm(Pfun2(x)-P*x) < 1e-12*norm(P*x), 'Incorrect Pfun2');
end

% Check modularity
k = sum(A,2);
Bx = A*x-k*k'*x/(2*m);
Bfun = mfunc_modularity(A);
Bfun2 = mfunc_modularity(@(x) A*x, n);
assert(norm(Bfun(x)-Bx) < 1e-12*norm(Bx), 'Incorrect Bfun');
assert(norm(Bfun2(x)-Bx) < 1e-12*norm(Bx), 'Incorrect Bfun2');

% Sanity check matrix scaling
[As,ab] = rescale_matrix(A);
x = randn(n,1);
y1 = As(x);
y2 = (A*x-ab(2)*x)/ab(1);
relerr = norm(y1-y2)/norm(y2);
assert(relerr < 1e-12, 'Inconsistent scaling behavior');
lambdas = (lambda-ab(2))/ab(1);
assert(all(abs(lambdas) < 1), 'Out-of-range scaled eigenvalues');

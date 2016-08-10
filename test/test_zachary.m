% Test basic functionality on the Zachary network
clear

% Load adjacency matrix and check basic
A = load_graph('zachary');
m = 78;  % Number of edges
n = 34;  % Number of nodes
assert(nnz(A) == 2*m,  'Wrong number of edges');
assert(size(A,1) == n, 'Wrong number of nodes');

% Compute full eigendecomposition for later sanity checks
[Q,D] = eig(full(A));
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

% Sanity check mfunc scaling
[As,ab] = rescale_mfunc(A);
x = randn(n,1);
y1 = As(x);
y2 = (A*x-ab(2)*x)/ab(1);
relerr = norm(y1-y2)/norm(y2);
assert(relerr < 1e-12, 'Inconsistent scaling behavior');
lambdas = (lambda-ab(2))/ab(1);
assert(all(abs(lambdas) < 1), 'Out-of-range scaled eigenvalues');

% Sanity check matrix scaling
[As,ab] = rescale_matrix(A);
y1 = As*x;
y2 = (A*x-ab(2)*x)/ab(1);
relerr = norm(y1-y2)/norm(y2);
assert(relerr < 1e-12, 'Inconsistent scaling behavior');
lambdas = (lambda-ab(2))/ab(1);
assert(all(abs(lambdas) < 1), 'Out-of-range scaled eigenvalues');

% Compare cheb moments for a given v
w = (Q'*x).^2;
c = moments_cheb(As, x, 10);
cref = zeros(10,1);
for k = 0:9
  cref(k+1) = w'*cos(k*acos(lambdas));
end
assert(norm(cref-c) < 1e-12, 'Inconsistent Cheb moments for fixed v');

% Compare cheb moments for a given v (second kind)
w = (Q'*x).^2;
c = moments_cheb(As, x, 10, 2);
cref = zeros(10,1);
for k = 0:9
 cref(k+1) = w'*(sin((k+1)*acos(lambdas))./sin(acos(lambdas)));
end
assert(norm(cref-c) < 1e-12, 'Inconsistent Cheb moments for fixed v');

% Compare cheb moments for DoS
[c,cs] = moments_cheb_dos(As, 1000, 10);
for k = 0:9
  cref(k+1) = sum(cos(k*acos(lambdas)));
end
relerr = abs(c-cref)./abs(c);
relerrb = cs./abs(c);
fprintf('--- Sanity check DoS moments (KPM) ---\n');
fprintf('%d: Relerr %e (vs %e for 95 pct CI): %d\n', ...
        [(0:9)', relerr, relerrb, 2*relerrb < relerr]');

[c,cs] = moments_lan_dos(As, 1000, 10);
relerr = abs(c-cref)./abs(c);
relerrb = cs./abs(c);
fprintf('--- Sanity check DoS moments (Lanczos) ---\n');
fprintf('%d: Relerr %e (vs %e for 95 pct CI): %d\n', ...
        [(0:9)', relerr, relerrb, 2*relerrb < relerr]');

% Compare cheb moments for LDoS
jnode = 12;
w = Q(jnode,:)'.^2;
[cx] = moments_cheb_ldosx(As, [jnode], 10);
[c,cs] = moments_cheb_ldos(As, 10);
c = c(:,jnode);
cs = cs(:,jnode);
for k = 0:9
  cref(k+1) = w'*cos(k*acos(lambdas));
end
relerr = abs(c-cref)./abs(c);
relerrx = abs(cx-cref)./abs(cx);
relerrb = cs./abs(c);
fprintf('--- Sanity check LDoS moments (KPM) ---\n');
fprintf('%d: Relerr %e (vs %e for 95 pct CI): %d\n', ...
        [(0:9)', relerr, relerrb, 2*relerrb < relerr]');
fprintf('--- Sanity check LDoS moments (vs purportedly exact) ---\n');
fprintf('%d: Relerr %e\n', ...
        [(0:9)', relerrx]');

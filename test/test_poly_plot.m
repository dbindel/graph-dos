% Load one of the Gleich examples
A = load_graph('gleich', 'musm-cc');

% Apply a filter
N = matrix_normalize(A);
c = filter_jackson(moments_delta(0.5, 500));
c(1) = c(1)/2;
pN = mfunc_cheb_poly(c,N);

% Pull out representative vectors (should revisit)
[X,~] = qr(randn(length(N), 100),0);
Y = pN(X);
[U,S,V] = svd(Y);
score = sum((U(:,1:2)*S(1:2,1:2)).^2, 2);

% Plot top group
[sscore,I] = sort(score, 'descend');
I = I(1:200);
As = A(I,I);
gplot_shatter(As,4);

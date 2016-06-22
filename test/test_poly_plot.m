% Load one of the Gleich examples
A = load_graph('gleich', 'musm-cc');

% Apply a filter
N = matrix_normalize(A);
c = filter_jackson(moments_delta(0.5, 500));
c(1) = c(1)/2;
pN = mfunc_cheb_poly(c,N);

% Pull out representative vectors (should revisit)
nprobe = 200;
thresh = 1e-3;
[V,D] = eig_rand1(randn(length(N),nprobe), pN, thresh);
[lambda,I] = sort(diag(D), 'descend');
V = V(:,I);

% Warn if it looks like we're missing something
if length(D) == nprobe
  warning('Did not converge to desired threshold\n');
end

% Sort out
score = sum(V(:,1:2).^2, 2);

% Plot top group
[sscore,I] = sort(score, 'descend');
I = I(1:100);
As = A(I,I);
gplot_shatter(As,4);

% K_5 with whiskers
A = [0 1 1 1 1 1 0 0;
     1 0 1 1 1 0 1 0;
     1 1 0 1 1 0 0 1;
     1 1 1 0 1 0 0 0;
     1 1 1 1 0 0 0 0;
     1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0];
 
A = matrix_normalize(A, 's');
[U,D,~] = svd(full(A));
d = diag(D);
hops = 4;
rtol = 1e-10;

sparse_vs = 0;
solved_vs = 0;

for j=1:length(A)
    lambda = d(j);
    v = U(:,j);
    nonzeros = find(abs(v) > rtol);
    if (length(nonzeros) > 3/4 * length(A))
        continue
    end
    center = nonzeros(1);
    sparse_vs = sparse_vs + 1;    
    [V,s] = recover_eigvec(A, lambda, center, hops, rtol);
    
    if (norm(A * V - lambda * V) < rtol && ~all(V))
        solved_vs = solved_vs + 1;
    end
end

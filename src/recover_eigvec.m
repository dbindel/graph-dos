% [V] = recover_eigenvec (A, lambda, center, max_hops)
%
% Determine the (sparse) eigenvector corresponding to lambda 
% that it is supported in neighborhood within h hops of specified center.
% 
% Inputs:
% - A:       given symmetric matrix in R^{n x n}
% - lambda:  eigenvalue of A
% - center:  a node known to support eigenvector, for center in [1, n]
% - h:       explore at most h hops from specified center

function [V, s] = recover_eigvec(A, lambda, center, h, rtol)
    n = length(A);
    
    % TODO how do we determine d?
    num_probe = 10;
    d = 10; % degree of function
    m = 2 * d; % number of moments
    
    [c,~] = moments_cheb_ldos (A, num_probe, m, 1);
    
    % Initialize visited array
    distance = zeros (1, n);
    distance(:) = intmax; distance(center) = 0;
    % Initialize queue for BFS
    queue = zeros (1, n);
    head = 1; last = 2;
    queue(head) = center;
    % Initialize support for eigenvector correponding to lambda
    s = zeros (n, 1);
    
    % Iterate through nodes in h-hop neighborhood
    while (head < last)
        node = queue(head);
        head = head + 1;
        d = distance(node);
        if (d == h)
            continue;
        end
        % Get the one-hop neighborhood
        neigh = A (node,:);
        for i=1:n
            if (neigh(i) ~= 0 && distance(i) == intmax || i == node)
                distance(i) = d + 1;
                m_seq = c(:,i);
                supp_eigs = freq_cheb(m_seq);
                if (any(abs(supp_eigs - lambda) < rtol))
                    queue(last) = i;
                    last = last + 1;
                    s(i) = 1;
                end
            end
        end
    end
    
    % Power iteration will converge to eigenvector v s.t. Av = lambda * v
    % and v and s have same support
    [~,V] = eigenbasis(A, lambda, s, rtol);

end
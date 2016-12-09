% [calc_eig] = freq_cheb (m_seq)
%
% Given Chebyshev moment sequence {m_seq} for LDoS on node k,
% return {calc_eig}, the eigenvalues with correponding
% eigenvectors supported on node k.

function [calc_eig] = freq_cheb (m_seq)

    m = length(m_seq); % number of moments
    d = m / 2; % degree

    H = hankel(m_seq(1:d+1), m_seq(d+1: m));
    [~,~,V] = svd(H);
    p = V(:, end);
    
    % Find the coefficients of the polynomial p given by
    % Chebyshev moments c_{km} and roots of p
    % p(x) = a_0 + a_1 x + \cdots + a_{d-1} x^{d-1} + a_d
    r = roots (flipud(p));
    % Fix roots so they are on the unit circle
    for j=1:length(r)
        n = norm([real(r(j)), imag(r(j))]); % should be approx 1
        r(j) = r(j) / n;
    end
    
    calc_eig = zeros(length(r), 1);
    
    for j=1:length(r)
        % r_j = e^{i w_j}
        % w_j = arccos (eig_j)
        w_j = log(r(j)) / sqrt(-1);
        eig_j = cos(w_j);
        assert (abs(imag(eig_j)) < 10^(-10));
        calc_eig(j,1) = real(eig_j);
    end
    calc_eig = abs(calc_eig);
end
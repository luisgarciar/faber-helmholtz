function sn = scalar_truncFS(n, M, N)
%SCALAR_TRUNCFS computes the truncated Faber series of f(z) = 1/z on a
%bratwurst set.
% 
% sn = SCALAR_TRUNCFS(n, M, N) takes as input the constants M, N of a bratwurst
% set (see [1]), and returns the coefficient vector sn of the n-th partial sum
% of the Faber series of 1/z on the bratwurst set.
% 
% References:
% [1] T. Koch and J. Liesen, The conformal ``bratwurst set'' and associated
%     Faber polynomials, Numer. Math., 86, (2000), pp. 173--191.
% 
% Author: Olivier Sète, Feb 14, 2016.
% 
% TO DO: If we need to speed things up: use representation of s_n(z) with Fhat's
% and use three-term recurrence.

alpha = - N - sqrt(N^2-1); % parameter for computing coeffs ak
ak = (M+alpha)/alpha;    % a_0
FF = get_Fpoly(n, M, N); % cell array of Faber polynomials F_0, ..., F_n
sn = ak * FF{1}; % s_0(z) = a_0 F_0(z)

% compute s_n(z) = sum_{k=0}^n a_k F_k(z)
for kk = 1:n
    ak = ak/alpha;
    sn = polyadd(sn, ak * FF{kk+1});
end

end
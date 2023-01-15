function s = FP1_pow(f, v, M, N)
%FP1_pow   Faber preconditioner of degree 1.
% 
%   s = FP1_pow(f, v, M, N) computes s1(A)*v, where s1(z) is the
%   truncated Faber series of 1/z of degree 1 on a `bratwurst' set.
%   f is a function handle with f(x) = A*x, v is a vector, and M, N are
%   constants for the `bratwurst' set.

% WARNING: Code below is specialised to n=1.

% Constants related to the `bratwurst' set:
alpha = - N - sqrt(N^2 - 1);
mu1 = 2*N;
% mu0 = 1;
nu1 = 2*(N - M);
nu0 = 2*(1 - M*N);

a1 = (M + alpha)/(alpha^2);
beta0 = a1 * (alpha - mu1 + nu0/nu1);
beta1 = a1 * nu1;

s = beta0*v;
s = s + beta1 * f(v);

s = f(s);   % A * s_1(A) * v

return

% Constant term of partial sum:
ak = (M + alpha)/alpha;     % a_0
s = (ak * (1 + nu0/(nu1*alpha) - mu1/alpha)) * v;

beta1 = nu1 * ak/alpha;     % coefficient of z
s = s + beta1 * f(v);          % s_1(A)*v

% s = f(s);   % A * s_1(A) * v

% F1 = nu1 * f(v) - mu1*v;    % Fhat_1(A) * v
% s = s + ak * F1;            % cst term + a_1 * Fhat_1(A) * v

end
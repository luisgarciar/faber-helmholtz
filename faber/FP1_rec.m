function s = FP1_rec(f, v, M, N)
%FP1_rec   Faber preconditioner of degree 1.
% 
%   s = FP1_rec(f, v, M, N) computes s = s1(A)*v, where s1(z) is the
%   truncated Faber series of 1/z of degree 1 on a `bratwurst' set.
%   f is a function handle with f(x) = A*x, v is a vector, and M, N are
%   constants for the `bratwurst' set.

% WARNING: Code below is specialised to n=1.

% Constants related to the `bratwurst' set:
alpha = - N - sqrt(N^2-1);
mu1 = 2*N;
% mu0 = 1;
nu1 = 2*(N-M);
nu0 = 2*(1-M*N);

% Constant term of partial sum:
quot = nu0/(nu1*alpha);
ak = (M + alpha)/alpha;     % a_0
s = (ak*(1 + quot)) * v;    % cst term of s1(A)*v

F1 = nu1 * f(v) - mu1*v;    % Fhat_1(A) * v
ak = ak/alpha;              % a_1
s = s + ak * F1;            % cst term + a_1 * Fhat_1(A) * v

s = f(s);                   % A * s1(A) * v

end
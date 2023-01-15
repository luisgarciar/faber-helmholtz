function s = FP2_pow(f, v, M, N)
%FP2_pow   Faber preconditioner of degree 2.
% 
%   s = FP2_pow(f, v, M, N) computes s2(A)*v, where s2(z) is the
%   truncated Faber series of 1/z of degree 2 on a `bratwurst' set.
%   f is a function handle with f(x) = A*x, v is a vector, and M, N are
%   constants for the `bratwurst' set.

% WARNING: Code below is specialised to n=2.

% Constants related to the `bratwurst' set:
alpha = - N - sqrt(N^2 - 1);
mu1 = 2*N;
% mu0 = 1;
nu1 = 2*(N - M);
nu0 = 2*(1 - M*N);

a2 = (M + alpha)/(alpha^3);
beta0 = a2*(alpha^2 - alpha*mu1 + alpha*nu0/nu1 + mu1^2 - 2 - (nu0/nu1)^2);
beta1 = a2*(alpha*nu1 + 2*nu0 - 2*mu1*nu1);
beta2 = a2*nu1^2;

s = beta0*v;
v = f(v);
s = s + beta1 * v;
v = f(v);
s = s + beta2 * v;

s = f(s);           % A * s_1(A) * v

end
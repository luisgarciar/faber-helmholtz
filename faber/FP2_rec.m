function s = FP2_rec(f, v, M, N)
%FP2_rec   Faber preconditioner of degree 1.
% 
%   s = FP2_rec(f, v, M, N) computes s = s2(A)*v, where s2(z) is the
%   truncated Faber series of 1/z of degree 2 on a `bratwurst' set.
%   f is a function handle with f(x) = A*x, v is a vector, and M, N are
%   constants for the `bratwurst' set.

% WARNING: Code below is specialised to n=2.


% get constants related to the bw set
alpha = - N - sqrt(N^2-1) ;
mu1 = 2*N ;
mu0 = 1 ;
nu1 = 2*(N-M) ;
nu0 = 2*(1-M*N) ;

% compute constant term of partial sum:
quot = nu0/(nu1*alpha);
ak = (M+alpha)/alpha;       % a_0
s = (ak*(2 - (1 + quot^3)/(1 + quot) )) * v; % cst term of sn(A)*v

% now we only need to add the terms a_k * Fhat_k(A)*v for k=1:n

F0 = 2*v ;                  % Fhat_0(A)*v
F1 = nu1 * f(v) - mu1*v;    % Fhat_1(A) * v
ak = ak/alpha;              % a_1
s = s + ak * F1;            % cst term + a_1 * Fhat_1(A) * v

F2 = f(nu1*F1 + nu0*F0) - mu1 * F1 - mu0 * F0;  % Compute Fhat_2
% rearranged, so that it needs only one matrix-vector multiplication
ak   = ak/alpha;            % a_2
s = s + ak * F2;            % shat_2 (+cst term)

s = f(s);                   % A * s_1(A) * v

end

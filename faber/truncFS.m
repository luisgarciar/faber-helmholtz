function snAv = truncFS(varargin)
% function snAv = TRUNCFS computes s_n(A)*v, where s_n(z) is a
% truncated Faber series of 1/z of degree n on a bw set.
%
% Three possibilities to call:
% 1) truncFS(A, v, n, M, N, 'mat')
% 2) truncFS(f, v, n, M, N, 'fun') where f is a function handle with 
%        f(x) = A*x.
% 3) truncFS(B, L, U, v, n, M, N, 'lu') where A = (LU)^-1 * B.
%
%
% Input:
%   A    -- square matrix (case 'mat')
%   f    -- function handle with f(v) = A*v (case 'fun')
%   A, L, U -- square matrices, where L, U are lu factors, and 
%   v    -- vector
%   n    -- degree at which Faber series of 1/z is truncated
%   M, N -- parameters of the bw set (get them from bw_map.m)
%
% Output:
%   snAv -- the value of s_n(A)*v
%
% For Details see [L. Garcia Ramos & O. Sete, Approximating the inverse on
% a bw set].
%

switch varargin{end}
    case 'mat'
        A = varargin{1} ;
        v = varargin{2} ;
        n = varargin{3} ;
        M = varargin{4} ;
        N = varargin{5} ;
        
        f = @(x) A*x ;
        
    case 'fun'
        f = varargin{1} ;
        v = varargin{2} ;
        n = varargin{3} ;
        M = varargin{4} ;
        N = varargin{5} ;
        
    case 'lu'
        A = varargin{1} ;
        L = varargin{2} ;
        U = varargin{3} ;
        v = varargin{4} ;
        n = varargin{5} ;
        M = varargin{6} ;
        N = varargin{7} ;
        
        f = @(x) U\(L\(A*x)) ;
        
    otherwise
        error('Last argument must specify mat, fun or lu.')
end

% get constants related to the bw set
alpha = - N - sqrt(N^2-1) ;
mu1 = 2*N ;
mu0 = 1 ;
nu1 = 2*(N-M) ;
nu0 = 2*(1-M*N) ;

% compute constant term of partial sum:
quot = nu0/(nu1*alpha) ;
ak = (M+alpha)/alpha ;     %% a_0
snAv = (ak*(2 - (1-(-quot)^(n+1))/(1+quot) )) * v ; %% cst term of sn(A)*v

if n == 0
    return
end

% now we only need to add the terms a_k * Fhat_k(A)*v for k=1:n

Fpp  = 2*v ;                  %% Fhat_0(A)*v
Fp   = nu1 * f(v) - mu1*v ;   %% Fhat_1(A) * v
ak   = ak/alpha ;             %% a_1
snAv = snAv + ak * Fp ;       %% cst term + a_1 * Fhat_1(A) * v
if n == 1
    return
end

Fcur = f(nu1*Fp + nu0*Fpp) - mu1 * Fp - mu0 * Fpp ;  %% Compute Fhat_2
% % rearranged, so that it needs only one matrix-vector multiplication
ak   = ak/alpha ;            %% a_2
snAv = snAv + ak * Fcur ;  %% shat_2 (+cst term)

for jj = 3:n
    Fpp  = Fp ;
    Fp   = Fcur ;
    Fcur = f(nu1*Fp + nu0*Fpp) - mu1*Fp - mu0*Fpp ;  %% F_jj
    ak   = ak / alpha ;
    snAv = snAv + ak * Fcur ;
end

end

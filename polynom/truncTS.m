function TPAv = truncTS(f,v,z0,n)
% function TPAv = truncTS(f,v,z0,n) computes T_n(A)*v, where T_n is the
% n-th order Taylor polynomial of 1/z at the point z0, and where f is a
% function handle with f(x) = A*x
% 
% Input:
%   f  -- function handle with f(x) = A*x
%   v  -- vector
%   z0 -- point for Taylor expansion
%   n  -- degree
% 
% Output:
%   TPAv -- the value of T_n(A)*v ;
% 
%   Author: Olivier Sete
%           TU Berlin, Institute of Mathematics
%   Version 0.1 - Nov 2015

% coefficients of the Taylor polynomial
quot = 1/z0 ;
ak = zeros(n+1,1) ;
for kk = 0:n
    ak(kk+1) = (-1)^kk * quot^(kk+1) ;  %% a_k = (-1)^k * (1/z0)^(k+1)
end

% Compute T_n(A)*v using Horner's rule
TPAv = ak(n)*v + ak(n+1)*f(v) - ak(n+1)*z0*v ; %% a_{n-1} + (z-z0)*a_n
for kk = (n-2):-1:0
    TPAv = ak(kk+1)*v + f(TPAv) - z0 * TPAv ;
end

end
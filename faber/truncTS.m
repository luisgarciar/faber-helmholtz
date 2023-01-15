function TnAv = truncTS(fA, v, n, z0)
%TRUNCTS    Compute Tn(A)*v, where Tn(z) is a Taylor polnomial of 1/z.
%
%   TnAv = TRUNCTS(fA, v, n, z0) computes Tn(A)*v, where fA is a function
%   handle implementing fA(x) = A*x, v is a vector, and Tn is the
%   Taylorpolynomial of 1/z of degree n and with center z0.

% Taylor expansion of f(z) = 1/z in z0:
%   f(z) = sum_{k=0}^\infty (-1)^k / z0^{k+1}   * (z-z0)^k

% Coefficient vector of T_n:
Tncoeff = 1/z0;
for kk = 1:n
    Tncoeff = polyadd(Tncoeff, (-1)^kk/z0^(kk+1)*poly(z0*ones(1,kk)));
end

% Compute Tn(A)*v:
curpow = v;
TnAv = Tncoeff(end)*curpow;
for kk=1:n
    curpow = fA(curpow);    % A^k = fA(A^(k-1))
    TnAv = TnAv + Tncoeff(end-kk)*curpow;
end

end
% 
% Test file -- investigate direct solve with Faber series for 1/z
% 
% We consider the CSL-preconditioned Helmholtz equation
%   A M^{-1} x = b.
% We write
%   x = (A M^{-1})^{-1} b = sum_{k=0}^Inf ak F_k(A)*b.
% The truncated series yields an approximation to the solution.  The converge of
% the series is slow and many terms are needed, already for small wavenumbers.
% We compare the result with GMRES, which is significantly faster.
% 
% Conclusion: 
% The Faber series of 1/z is not suitable to solve Ax=b directly in the
% Helmholtz setting.


k  = 40;        % Wavenumber
npc = 1;        % Number of interior points in coarsest grid in one dim.
bc = 'dir';     % Boundary condition.
ppw = 15;       % Number of points per wavelength.
eps = 0.5*k^2;  % Shift for CSL.

% Number of points in finest grid and number of levels:
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);

% Matrices (1D):
A = helmholtz(k,0,npf,bc);
M = helmholtz(k,eps,npf,bc);

% LU-decomposition of M:
[L, U] = lu(M);

%Preconditioned matrix (with DBC)
f = @(x) A*(U\(L\x));   % f(x) = A * M^{-1} * x

% Parameters for the bratwurst shaped set
lambda = -1;
phi = 0.1 * pi;
eps_thick = 0.005;
[~, ~, ~, Mbw, Nbw] = bw_map(lambda, phi, eps_thick);

% Construct linear system:
x_exact = ones(size(A, 2),1);   % Solution and
b = f(x_exact);                 % right hand side.

% Note:  To speed up the computation, the partial sums should be computed step
% by step (just adding the next factor), not recomputed for each degree.

% Compute solution by partial sumes:
for deg = 1:600
    x = truncFS(f, b, deg, Mbw, Nbw, 'fun');    % Compute s_deg(A M^{-1})*b
    err(deg) = norm(x-x_exact);
    res(deg) = norm(b - f(x));
end

% plot parameters
FS = 22;

% plot the error
figure(1)
semilogy(err)
largefiglabels
title('Error |x-x_{exact}|')
xlabel('degree of Faber series')
ylim([10^-14  10^2])




% Solve with GMRES
restart = [];
tol = 10^-12;
maxit = size(A, 1);
[x_gmres, ~, ~, iter, resvec] = gmres(f, b, restart, tol, maxit);

err_gmres = norm(x - x_gmres)

% Plot relative residuals
figure(2)
ph_fs = semilogy(res ./ res(1));
hold on
ph_gmres = semilogy(resvec ./ resvec(1), 'r-');
hold off
largefiglabels
title('Relative residual norms')
xlabel('deg of Faber S / number of iteration in GMRES')
ylim([10^-14  10^2])
legend('rel res (Faber series)', 'rel res (GMRES)')

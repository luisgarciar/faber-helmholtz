 
% Test file : compare number of iterations and computation time
% for solving the 1-d Helmholtz equation preconditioned with a shifted
% Laplacian.  We compare solving with GMRES and BiCGstab without further
% pre-conditioning, and with "preconditioning" using a truncated Faber
% series for 1/z.
% (Using the new implementation.)
% 

k  = 400;    %wavenumber
np = ceil( 12 * k / pi) ;
b  = 0.5;  %complex shift 
bc = 'dir' ; %type of boundary conditions

%1-D Helmholtz
A_d = helmholtz(k,np,bc);

%1-D Shifted Laplacian
M_d = shift_laplace(k,1,b,np,bc);

%Preconditioned matrix (with DBC)
S_d = M_d\A_d;
f = @(x) M_d \ (A_d*x) ;

% Parameters for the bratwurst shaped set
lambda = -1 ;
phi = 0.1 * pi ;
eps_thick = 0.005 ;
[psi, ~, capacity, M, N] = bw_map(lambda, phi, eps_thick) ;

b = ones(np, 1) ;

% parameters for GMRES and BiCGstab
restart = [] ;  %% number of iter before restart, [] = no restart.
tol = 1e-12 ;   %% tolerance for relative residual

tic
[X_gmres,~,~,iter_gmres,resvec_gmres] = gmres(f, b, restart, tol, np) ;
time_gmres = toc

tic
[X_bicg,~,~,iter_bicg,resvec_bicg] = bicgstab(S_d, b, tol, np) ;
     % with f instead of S_d the method diverges.
time_bicgstab = toc
% note:
% bicg stores in resvec every iteration + half iteration.
% Hence the adjustement in the plots below

deg = 5 ;  %% degree where to truncate the Faber series of 1/z

tic
snb = truncFS(f,b,deg, M,N,'fun') ;  %% compute new rhs: s_n(S_d) * b
[X_FS,~,~,iter_FS,resvec_FS] = gmres(@(x) truncFS(f,f(x),deg,M,N,'fun'),...
    snb, restart, tol, np) ;
time_FS = toc

tic
snb = truncFS(f, b, deg, M, N, 'fun') ;
[X_bicg_FS,~,~,iter_bicg_FS, resvec_bicg_FS] = ...
    bicgstab(@(x) truncFS(f,f(x),deg,M,N, 'fun'), snb, tol, np) ;
time_bicgstab_FS = toc


figure(1)
semilogy(resvec_gmres/resvec_gmres(1), 'k-')
hold on
semilogy(resvec_FS/resvec_FS(1), 'k--')
semilogy( resvec_bicg(1:2:end) / resvec_bicg(1), 'r-')
semilogy( resvec_bicg_FS(1:2:end) / resvec_bicg_FS(1), 'r--')
hold off


%Experiments with Faber acceleration and BiCGstab

% Experiments with the CSL-preconditioned 1D Helmholtz equation.
% We compare timing and the iteration numbers for solving
%   AM^{-1} x = b
% with gmres, for different implementations of M^{-1}Ax:
%   - ilu: U\(L\(Ax))
%   - multi-grid (TO DO)
% 
% Further, we do the same when using polynomial acceleration with a truncated
% Faber series s_n(z) of 1/z on a bw set. Then the system to solve is
%   s_n(M^{-1}A) M^{-1}A x = s_n(M^{-1}A) b.
% 
% OBSERVATIONS: (Dirichlet boundary conditions)
% - Generally ilu is faster than the two other methods (with or without poly
%   acceleration).
% - With truncated Faber series less iterations are needed.
% - With truncated Faber series less time is needed starting at k = 50.

% Setup list of wavenumbers
clear all
close all

%kk  = [10 20 30 40 50]'
%npp = [31 63 127 127 255]'
pmin = 20;       %min points per wavelenght
k    = 250;   %% run this when testing changes in the code
lev = ceil(log2(k*pmin/2*pi)); npf = 2^lev-1; %number of 1D interior gridpoints

b1=1; b2=0.5;  %complex shift
bc   = 'dir' ;   %type of boundary conditions
dim  = 1;     
A    = helmholtz(k,npf,bc);
M    = shift_laplace(k,b1,b2,npf,bc);
b = zeros(length(A),1); b(ceil(length(A)/2),1)=1;
z = zeros(length(A),1);

%Parameters for BiCG-STAB
tol    = 1e-4;
maxit = length(A);

% Multigrid Parameters
npre = 2; npos = 2; w  = 2/3; smo = 'wjac'; numcycles=1;
       

% Creating the matrices and intergrid operators on all levels
 [SLgrid_matrices,SLgrid_smooth,restrict,interp] = mgsm_setup(M,lev,bc,dim);

tic
[iL, iU] = ilu(M); % iL * iU approx M
time_ilu1 = toc;


%Setup iLU and MG preconditioners

% Setting the SL preconditioner Minv
Minv = @(x)feval(@Vcycle2,SLgrid_matrices,SLgrid_smooth,...
                 restrict,interp,z,x,npre,npos,w,smo,numcycles); 

%Setup iLU and MG preconditioners
f_ilu = @(x) A*(iU\(iL\x));
f_mg  = @(x) A*feval(Minv,x); 

%  == Solve AM^{-1} x = b    
%  with f_ilu
tic
[~, ~, ~, iter_ilu, resvec_ilu] = bicgstab(f_ilu, b, tol, maxit);
time_ilu = toc;

% with f_mg~
tic
[~, ~, ~, iter_mg, resvec_mg] = bicgstab(f_mg, b, tol, maxit);
time_mg = toc;
         
%%  Polynomial acceleration with truncated Faber series s_deg then is
% s_deg( M^{-1}A ) M^{-1}A x = s_deg( M^{-1}A )b
% Compare with same computations as before.
deg = 2; % degree of truncated Faber series of 1/z.
     
%Parameters for the bratwurst shaped set
lambda    = -1; % so that 0 is not in the bw set.
phi       = pi/2;  % or phi = 0.1*pi; ??
eps_thick = 0.005;
[~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);

%Setup iLU_FS and MG_FS preconditioners
f_ilu_FS = @(x)truncFS(f_ilu,f_ilu(x),deg,M_bw,N_bw,'fun');
f_mg_FS  = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');

% with f_ilu_FS
tic
[~,~,~,iter_ilu_FS,resvec_ilu_FS] = bicgstab(f_ilu_FS,b,tol,maxit);
time_ilu_FS = toc;
   
% with f_mg
tic
[~,~,~,iter_mg_FS,resvec_mg_FS] = bicgstab(f_mg_FS,b,tol,maxit);
time_mg_FS = toc;

%% % % Plot time
 
% figure(1)
% plot(wavenum, time_lu, 'k--')
% hold on
% %plot(wavenum, time_ilu,'g-')
% plot(wavenum, time_mg, 'k-')
% plot(wavenum, time_lu_FS, 'b--')
% %plot(wavenum, time_ilu_FS, 'r-')
% plot(wavenum, time_mg_FS, 'b-')
% hold off
% ylabel('time (s)')
% xlabel('wavenumber k')
% legend('f_{lu}', 'f_{ilu}', ...
%      'f_{lu} and poly acc with FS', ...
%      'f_{ilu} and poly acc with FS', 'Location','NorthWest')
% title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),') and deg=',...
%     num2str(deg)])
% FS = 22 ; % font size
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'FontSize',FS);

%break;

%% Plot iteration numbers
% figure(2)
% plot(wavenum, iter_lu(:,2), 'k--')
% hold on
% plot(wavenum, iter_ilu(:,2), 'g-')
% plot(wavenum, iter_mg(:,2),    'k-')
% plot(wavenum, iter_lu_FS(:,2), 'b--')
% plot(wavenum, iter_ilu_FS(:,2), 'r-')
% plot(wavenum, iter_mg_FS(:,2), 'b-')
% hold off
% ylabel('gmres iterations')
% xlabel('wavenumber k')
% legend('f_{lu}', 'f_{ilu}','f_mg', ...
%     'f_{lu} and poly acc with FS', ...
%     'f_{ilu} and poly acc with FS', 'Location','NorthWest');
% title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%     num2str(deg)])
% % FS = 22 ; % font size
% % set(gca,'LooseInset',get(gca,'TightInset'))
% % set(gca,'FontSize',FS);
% 
% % figure(1)
% % print('-depsc2', ['Helm1D_',bc,'_FS_timing_in_k.eps'])
% % figure(2)
% % print('-depsc2', ['Helm1D_',bc,'_FS_iterations_in_k.eps'])
% 
% break

iter_ilu
time_ilu 
iter_ilu_FS
time_ilu_FS 

iter_mg
time_mg
iter_mg_FS
time_mg_FS   

% Plot the relative residual curves FOR ONE WAVENUMBER k
figure(1)
semilogy(0:iter_ilu(1),resvec_ilu(1:2:end)/resvec_ilu(1), 'k-.')
hold on
semilogy(0:iter_mg(1),resvec_mg(1:2:end)/resvec_mg(1), 'b-.')
semilogy(0:iter_ilu_FS(1),resvec_ilu_FS(1:2:end)/resvec_ilu_FS(1), 'k-+')
semilogy(0:iter_mg_FS(1), resvec_mg_FS(1:2:end)/resvec_mg_FS(1), 'b-+')
hold off
ylabel('relative residual')
xlabel('iteration')
legend('f_{ilu}','f_mg', ...
    'f_{ilu} and poly acc with FS','f_{mg} and poly acc with FS')
title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),')'])
FS = 22 ; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS);

% print('-depsc2', 'Helm1D_rr_FS_different_MinvA.eps')



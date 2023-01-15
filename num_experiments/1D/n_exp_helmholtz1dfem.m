%% Experiments with the CSL-preconditioned 1D Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    A M^{-1}x = b
% with gmres, for different implementations of M^{-1}Ax:
%   - "exact": M\(Ax)
%   - lu: U\(L\(Ax))
%   - multi-grid
% 
% Further, we do the same when using polynomial acceleration with a truncated
% Faber series s_n(z) of 1/z on a bw set. Then the system to solve is
%   s_n(M^{-1}A) M^{-1}A x = s_n(M^{-1}A) b.
% 

% Setup list of wavenumbers
%wavenum = [10:10:90, 100:100:500];
%wavenum  = [100:100:300];
%wavenum = [10:10:40];
%wavenum = [200 300];
wavenum = 200;       %% run this when testing changes in the code
ppw = 12;            %min points per wavelenght
bc = 'som';
npc = 3;              %number of points coarsest grid
eps = 0.5*(wavenum.^2);     %complex shift
dim  = 1;     
npre = 2; npos = 2; numit = 1; %number of pre, postsmoothing steps (multigrid)
smo= 'wjac'; w  = 2/3;         %damping parameter for wJacobi 

%Parameters for GMRES
restart = [];% [];
tol     = 1e-8;
maxit   = 300;

%memory allocation for time, residuals
time_lu     = zeros(length(wavenum),1);
time_lu_FS  = zeros(length(wavenum),1);
time_mg     = zeros(length(wavenum),1);
time_mg_FS  = zeros(length(wavenum),1);
time_mgl    = zeros(length(wavenum),1);
time_mg_FSl = zeros(length(wavenum),1);


iter_lu     = zeros(length(wavenum),2);
iter_lu_FS  = zeros(length(wavenum),2);
iter_mg     = zeros(length(wavenum),2);
iter_mg_FS  = zeros(length(wavenum),2);
iter_mg_l   = zeros(length(wavenum),2);
iter_mg_FSl = zeros(length(wavenum),2);

%% Iteration phase
for j = 1:length(wavenum)
    k = wavenum(j); %wavenumber
    npc = 3;
    [npf,numlev] = fem_npc_to_npf(npc,k,ppw); 
    A = helmholtzfem(k, npf, 0,bc);      %1D Helmholtz matrix
    M = helmholtzfem(k, npf, eps(j),bc); %CSL with complex shift eps
    
    op_type = 'gal';
     
    [SLgrid_matrices,SLgrid_split,restrict,interp] = ...
         mg_setupfem(k,eps(j),op_type,npc,numlev,dim,bc); %multigrid setup for shifted Laplacian
    z = zeros(length(M),1);
   
    Minv_lu   = @(x) U\(L\x);
    AMinv_lu  = @(x) A*feval(Minv_lu,x);
    Minv_mg   = @(x) feval(@Vcyclefem,SLgrid_matrices,SLgrid_split,restrict,interp,z,x,npre,npos,w,smo,numit); 
    AMinv_mg  = @(x) A*feval(Minv_mg,x);
    
    h = 1/npf; x = h*(1:1:npf)';
    b = ones(npf,1); b(npf) = 0.5; b = h*b;    %right hand side
    ex = exact_sol(k,x);
    
    % %== Ax=b is TOO SLOW, already for k=60: 1s, k=100: 3s on my Laptop. We may run
    % % this *once* if we must, but not for usual experiments.
%      tic
%      [x_unpr, ~, ~, iter_unprec, resvec_unprec] = gmres(A, b, restart, tol, npf);
%      time_unprec = toc;
    
    % == Solve AMinv x = b
    % with AMinv_lu
%     tic
%     [x_lu, ~, ~, iter_lu(j,:), resvec_lu] = gmres(AMinv_lu, b, restart, tol, maxit);
%     time_lu(j,1) = toc;
%     

    % with AMinv_mg
     tic
     [x_mg, ~, ~, iter_mg(j,:), resvec_mg] = gmres(AMinv_mg, b, restart, tol, maxit);
     time_mg(j,1) = toc;
     
      % with Minv_mgA
      tic
      [x_mgl, ~, ~, iter_mgl(j,:), resvec_mgl] = gmres(A, b, restart, tol, maxit,Minv_mg);
      time_mgl(j,1) = toc;
%     
     
     %  Polynomial acceleration with truncated Faber series s_deg then is
     % s_deg( M^{-1}A ) M^{-1}A x = s_deg( M^{-1}A )b
     % Compare with same computations as before.
     
     deg = 1; % degree of truncated Faber series of 1/z.
     
     %Parameters for the bratwurst shaped set
     lambda    = -1;    % so that 0 is not in the bw set.
     phi       = pi/3;  % opening angle
     eps_thick = 0.1;   %0.005;
     [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
     
      %Setting the Polynomial preconditioners: Faber with LU, Faber with MG
      AMinv_mg_FS  = @(x) AMinv_mg(truncFS(AMinv_mg,x,deg,M_bw,N_bw,'fun'));
      %AMinv_lu_FS  = @(x) AMinv_lu(truncFS(AMinv_lu,x,deg,M_bw,N_bw,'fun'));
      
      % with mg+FS
      tic
      [x_mg_FS, ~, ~, iter_mg_FS(j,:),resvec_mg_FS] = gmres(AMinv_mg_FS, b, restart, tol, maxit);
      time_mg_FS(j,1) = toc;
      
%       tic
%       [x_lu_FS, ~, ~, iter_lu_FS(j,:),resvec_lu_FS] = gmres(AMinv_lu_FS, b, restart, tol, maxit);
%       time_lu_FS(j,1) = toc;

 end % of going through different wavenumbers

 %time_unprec
 %time_lu
 %time_lu_FS 
 time_mgl
 time_mg
 time_mg_FS
 
 %iter_unprec
%  iter_lu
%  iter_lu_FS
 iter_mgl
 iter_mg
 iter_mg_FS
 
 gain_time = 100*(time_mg-time_mg_FS)/time_mg 
 gain_iter = 100*(iter_mg-iter_mg_FS)/iter_mg
 
 h         = 1/npf; grid = h*(1:1:npf)';
 u_ex      = exact_sol(k,x);
 u         = A\b;  plot(grid,real(u),'k-',grid,real(u_ex),'r-');
 relerr    = norm(u-u_ex)/norm(u_ex)
 mglrelerr = norm(x_mgl-u_ex)/norm(u_ex)

%% % %% Plot time

% figure(1)
% plot(wavenum, time_mg, 'k-')
% hold on
% plot(wavenum, time_mg_FS, 'b-')
% hold off
% ylabel('time (s)')
% xlabel('wavenumber k')
% legend('MG', 'MG+FS', 'Location','NorthWest');
%title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),') and deg=',...
   % num2str(deg)])
%FS = 22 ; % font size
%set(gca,'LooseInset',get(gca,'TightInset'))
%set(gca,'FontSize',FS);

% break;
% 
% 
% % %% Plot iteration numbers
% % figure(2)
% % plot(wavenum, iter_lu(:,2), 'k--')
% % hold on
% % plot(wavenum, iter_ilu(:,2), 'g-')
% % plot(wavenum, iter_mg(:,2),    'k-')
% % plot(wavenum, iter_lu_FS(:,2), 'b--')
% % plot(wavenum, iter_ilu_FS(:,2), 'r-')
% % plot(wavenum, iter_mg_FS(:,2), 'b-')
% % hold off
% % ylabel('gmres iterations')
% % xlabel('wavenumber k')
% % legend('f_{lu}', 'f_{ilu}','f_mg', ...
% %     'f_{lu} and poly acc with FS', ...
% %     'f_{ilu} and poly acc with FS', 'Location','NorthWest');
% % title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
% %     num2str(deg)])
% % % FS = 22 ; % font size
% % % set(gca,'LooseInset',get(gca,'TightInset'))
% % % set(gca,'FontSize',FS);
% % 
% % % figure(1)
% % % print('-depsc2', ['Helm1D_',bc,'_FS_timing_in_k.eps'])
% % % figure(2)
% % % print('-depsc2', ['Helm1D_',bc,'_FS_iterations_in_k.eps'])
% % 
% % break
% % 
% % % Plot the relative residual curves FOR ONE WAVENUMBER k
% % figure(1)
% % semilogy(0:iter_lu(2), resvec_lu/resvec_lu(1), 'k-')
% % hold on
% % semilogy(0:iter_ilu(2), resvec_ilu/resvec_ilu(1), 'k-')
% % semilogy(0:iter_mg(2), resvec_mg/resvec_mg(1), 'k-')
% % semilogy(0:iter_ex_FS(2), resvec_ex_FS/resvec_ex_FS(1), 'b-')
% % semilogy(0:iter_lu_FS(2), resvec_lu_FS/resvec_lu_FS(1), 'b-')
% % semilogy(0:iter_ilu_FS(2), resvec_ilu_FS/resvec_ilu_FS(1), 'b-')
% % semilogy(0:iter_mg_FS(2), resvec_mg_FS/resvec_mg_FS(1), 'b-')
% % hold off
% % ylabel('relative residual')
% % xlabel('iteration')
% % legend('f_{ex} = M^{-1}A*x (plain gmres)', 'f_{lu}', 'f_{ilu}', ...
% %     'f_{ex} and poly acc with FS', 'f_{lu} and poly acc with FS', ...
% %     'f_{ilu} and poly acc with FS')
% % title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),')'])
% % % FS = 22 ; % font size
% % % set(gca,'LooseInset',get(gca,'TightInset'))
% % % set(gca,'FontSize',FS);
% % 
% % % print('-depsc2', 'Helm1D_rr_FS_different_MinvA.eps')
% Experiments with the CSL-preconditioned 2D Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    AM^{-1}x = b
% with gmres, for different implementations of A*M^{-1}*x:
%   - multi-grid
%
% The 2D Sommerfeld problem is discretized with finite differences.
%
% Further, we do the same when applying the Faber preconditioner s_n:
%   AM^{-1} s_n(AM^{-1}) x = b.

% Computation: this file
% Postprocessing: mp3_postprocess_iter_time

close all
clear vars;

%% Setup parameters
%Setup list of wavenumbers
wavenum = [5 10 20]; %% warmup
% wavenum = 50:50:400;
degree  = 1:2;

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 200;  %maximum number of gmres cycles
numruns = 1;    %number of runs for every experiment (to average later)

num_dofs    = zeros(length(wavenum), 1);
time_mg     = zeros(length(wavenum), 1);
time_mg_FS  = zeros(length(wavenum), length(degree));
iter_mg     = zeros(length(wavenum), 2);
iter_mg_FS  = zeros(length(wavenum), length(degree), 2);

for s = 1:numruns
    
    for kk = 1:length(wavenum)
        %wavenumber
        k = wavenum(kk);
        factoreps = 0.5;
        poweps    = 2;
        eps       = factoreps*k^poweps;
        npc       = 4;   %number of points in the coarsest grid in 1D
        ppw       = 12;  %pollution free grid
        bc        = 'som1';
        %number of points in finest grid and number of levels
        [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
        
        %2D
        dim  = 2;
        
        %Helmholtz matrix in 2D
        A = helmholtz2(k,0,npf,npf,bc);
        
        %Multigrid structure for shifted Laplacian
        op_type = 'gal';
        [mg_mat_SL,mg_split,restr,interp]= mg_setup(k,eps,op_type,npc,numlev,bc,dim);
        
        %Shifted Laplace matrix
        Aeps = mg_mat_SL{1};
        
        %Right hand side
        b  = zeros(length(A),1); ind = floor(length(b)/2);  b(ind)=1;
        x0 = zeros(size(b));
        
        %multigrid and gmres parameters
        npre  = 1; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 1;
        f_mg  = @(v) A*feval(@Fcycle,mg_mat_SL,mg_split,restr,interp,x0,v,npre,npos,w,smo,1);
        
        %run gmres only for shifted Laplacian
        tic
        [~,flag_mg,~,iter_mg(kk,:),resvec_mg] = gmres(f_mg, b, restart, tol, maxit);
        time_mg(kk) = time_mg(kk) + toc;
        
        for deg = 1:length(degree)
            %Setup of Polynomial acceleration with truncated Faber series
            %Parameters for the bratwurst shaped set
            lambda    = -1; % so that 0 is not in the bw set.
            phi       = pi/2;  % or phi = 0.1*pi; ??
            eps_thick = 0.005;
            [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
            
            %Setup  MG_FS preconditioners
            f_mg_FS = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');
            
            % with f_mg_FS
            tic
            [~,flag_mg_FS,~,iter_mg_FS(kk,deg,:),resvec_mg_FS] = gmres(f_mg_FS,b,restart,tol,maxit);
            time_mg_FS(kk,deg) = time_mg_FS(kk,deg) + toc;
            
        end %of going through different degrees
    end %of going through different wavenumbers
end %of going through different runs


if (isempty(restart))
    restart = 0;
end

totiter_mg = (iter_mg(:,1)-1)*restart + iter_mg(:,2);
mvops_mg   = totiter_mg;

totiter_mg_FS = (iter_mg_FS(:,:,1)-1)*restart + iter_mg_FS(:,:,2);
mvops_mg_FS   = totiter_mg_FS * diag(degree+1);

time_mg_FS = time_mg_FS/numruns;
time_mg    = time_mg/numruns;

save('mp3_results_cluster', 'wavenum', 'dim', 'degree', 'restart', 'bc',...
    'time_mg', 'time_mg_FS', 'totiter_mg', 'totiter_mg_FS', ...
    'mvops_mg', 'mvops_mg_FS')

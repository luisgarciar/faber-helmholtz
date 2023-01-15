% Experiments with the CSL-preconditioned 2D Helmholtz equation.
% We compare timing and the iteration numbers for solving
%    AM^{-1}x = b
% with gmres, for different implementations of A*M^{-1}*x:
%   - multi-grid
%
% The 2D Sommerfeld problem is discretized with finite differences
%
% Further, we do the same when using polynomial acceleration with a truncated
% Faber series s_n(z) of 1/z on a bw set. Then the system to solve is
%   AM^{-1}s_n(AM^{-1}) x = b.
%
% OBSERVATIONS: (Sommerfeld boundary conditions)

close all
clear all;
save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.

%% Setup parameters
%Setup list of wavenumbers
%wavenum  = [20:20:120];

%Parameters of shifted Laplacian
k = 80;
factoreps = 0.6;
poweps    = 2;
eps       = factoreps*k^poweps;
npc       = 1;    %number of interior points in the coarsest grid in 1D
ppw       = 12; %pollution free grid
bc        = 'som1';

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 100; % maximum number of gmres cycles
numruns = 1; %number of runs for every experiment (to average later)

%number of points in finest grid and number of levels
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);

%2D
dim = 2;

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
[~,flag_mg,~,iter_mg,resvec_mg] = gmres(f_mg, b, restart, tol, maxit);

if (isempty(restart))
    rest = 0;
else
    rest = restart;
end
totiter_mg = (iter_mg(1)-1)*rest + iter_mg(2);
mvops_mg   = 0:totiter_mg;

figure(1)
semilogy(mvops_mg, resvec_mg, 'k-','Linewidth',2)
hold on
ylabel('Relative residual')
xlabel('Number of matrix-vector products')
leg{1} = 'CSL(MG)';


FS = 16; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)

degree = 1:3;
linecolor = {'b','r','m'};

%% Run Experiment
for deg = 1:length(degree)
    
    %Setup of Polynomial acceleration with truncated Faber series
    %Parameters for the bratwurst shaped set
    lambda    = -1; % so that 0 is not in the bw set.
    phi       = pi/2;  % or phi = 0.1*pi; ??
    eps_thick = 0.005;
    [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
    
    %Setup MG_FS preconditioner
    f_mg_FS  = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');
    
    for t=1:numruns
        
        % with f_mg_FS
        tic
        [X_mg_FS,flag_mg_FS,~,iter_mg_FS,resvec_mg_FS] = gmres(f_mg_FS,b,restart,tol,maxit);
        
        if (isempty(restart))
            rest = 0;
        else
            rest = restart;
        end
        
        % Count matrix-vector operations
        % For GMRES on A*x = b, each step requires one multiplication by the matrix A.
        % For GMRES on p(A)A*x = b, each step requires one multiplication by the matrix
        % p(A)A, i.e., 1 + deg(p) multiplications by A.
        % Applying this to the matrix M^{-1} A in our problem, we see that MG
        % has ITER multiplications by M^{-1}A, i.e., M^{-1} is applied ITER times.
        % For MG+FS, we have ITER * (1+d) multiplications by M^{-1}A and
        % applications of M^{-1}.
               
        totiter_mg_FS = (iter_mg_FS(1)-1)*rest + iter_mg_FS(2);
        mvops_mg_FS   = (deg+1)*(0:totiter_mg_FS);
        
        semilogy(mvops_mg_FS, resvec_mg_FS, 'color',linecolor{deg},...
            'linestyle','-','Linewidth',2)
        hold on
        
        leg{deg+1} = strcat('FP(',num2str(deg),')+CSL(MG)');
        
        legend(leg,'Location','NorthEast')
    end
    
end % of going through different wavenumbers


%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the matlab path
wn     = num2str(k);
dimn   = num2str(dim);
degmin = num2str(min(degree));
degmax = num2str(max(degree));
restt  = num2str(rest); 
bdc    = num2str(bc);

plot_mvops_tex = strcat('mvops_vs_deg_k=',wn,'_',dimn,'D_',...
    bdc,'deg',degmin,'to',degmax,'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');

%generate file name, for saving in tex_files/figures/new_exp
path          = fullfile(currentpath,'..','..','..','..','tex_files',...
    'figures','new_exp');

plot_file_tex = fullfile(path,plot_mvops_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    plot_time_eps = strcat('mvops_vs_deg_k=',wn,'_',dimn,'D_',...
                            bdc,'deg',degmin,'to',degmax,...
                            'restart',restt,'.eps');
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end


profile off

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

% Computation and postprocessing.

close all
clear vars;
save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.
LW = 'LineWidth'; lw = 1;

%% Setup parameters
%Setup list of wavenumbers
wavenum = [5 10 20 40 60 80 100 120 140];
wavenum = [5, 20, 50, 80, 110, 140];
degree  = 1:2;

%Parameters for GMRES
restart = [];
% restart = 20;
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


%% Plot timings and save the plot in tikz (.tex) and .eps formats
figure(1)
plot(wavenum, time_mg, 'k-o', LW, lw);
hold on
plot(wavenum, time_mg_FS(:,1), 'b--x', LW, lw);
plot(wavenum, time_mg_FS(:,2), 'r-.+', LW, lw);
hold off

ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL(MG)', ...
    'FP(1)+CSL(MG)',...
    'FP(2)+CSL(MG)',...
    'Location','NorthWest','FontSize',16,'interpreter','latex');


%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the matlab path
kmin    = num2str(min(wavenum));
kmax    = num2str(max(wavenum));
dimn    = num2str(dim);
degmin  = num2str(min(degree));
degmax  = num2str(max(degree));
restt   = num2str(restart);
bdc     = num2str(bc);

plot_time_tex = strcat('time_vs_k_',kmin,'to',kmax,...
    '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
    'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');

%generate file name, for saving in tex_files/figures/new_exp
% path          = fullfile(currentpath,'..','..','..','..','tex_files',...
%     'figures','new_exp');
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
path = fullfile(currentpath, 'tex_files', 'figures', 'new_exp');


plot_file_tex = fullfile(path,plot_time_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex, 'width', '\figWidth', ...
        'extraaxisoptions',...
        ['xlabel style={font=\scriptsize, at={(axis description cs:0.5,-0.05)},anchor=north},', ...
        'ylabel style={font=\scriptsize, at={(axis description cs:-0.05,.5)},anchor=south},', ...
        'legend style={font=\scriptsize},']);
    
    %Save the plot as .eps
    plot_time_eps = strcat('time_vs_k_',kmin,'to',kmax,...
        '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
        'restart',restt,'.eps');
    
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end


%% Plot iteration numbers
figure(2)
plot(wavenum, totiter_mg, 'k-o', LW, lw);
hold on
plot(wavenum, totiter_mg_FS(:,1), 'b--x', LW, lw);
plot(wavenum, totiter_mg_FS(:,2), 'r-.+', LW, lw);
hold off

ylabel('Number of iterations')
xlabel('Wavenumber')
legend('CSL(MG)', ...
    'FP(1)+CSL(MG)',...
    'FP(2)+CSL(MG)',...
    'Location','NorthWest','FontSize',16,'interpreter','latex')


%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the matlab path
kmin    = num2str(min(wavenum));
kmax    = num2str(max(wavenum));
dimn    = num2str(dim);
degmin  = num2str(min(degree));
degmax  = num2str(max(degree));
restt   = num2str(restart);
bdc     = num2str(bc);

plot_iter_tex = strcat('iter_vs_k_',kmin,'to',kmax,...
    '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
    'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');

%generate file name, for saving in tex_files/figures/new_exp
% path          = fullfile(currentpath,'..','..','..','..','tex_files',...
%     'figures','new_exp');
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
path = fullfile(currentpath, 'tex_files', 'figures', 'new_exp');


plot_file_tex = fullfile(path,plot_iter_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex, 'width', '\figWidth', ...
        'extraaxisoptions',...
        ['xlabel style={font=\scriptsize, at={(axis description cs:0.5,-0.05)},anchor=north},', ...
        'ylabel style={font=\scriptsize, at={(axis description cs:-0.09,.5)},anchor=south},', ...
        'legend style={font=\scriptsize},']);
    
    %Save the plot as .eps
    plot_iter_eps = strcat('iter_vs_k_',kmin,'to',kmax,...
        '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
        'restart',restt,'.eps');
    
    plot_file_eps = fullfile(path,plot_iter_eps);
    print('-depsc2',plot_file_eps)
end


%% Export .tex table with the timings for several preconditioning methods
%Table with comparison data of MG
table1_data = [wavenum.',totiter_mg,mvops_mg,time_mg,...
               totiter_mg_FS(:,1),mvops_mg_FS(:,1),time_mg_FS(:,1),...
               totiter_mg_FS(:,2),mvops_mg_FS(:,2),time_mg_FS(:,2)];

numrows      = size(table1_data,1);
numcols      = size(table1_data,2);
if ( restart == 0 )
    tableCaption = strcat('Model problem 3 with full GMRES: Number of GMRES iterations, matrix-vector products and timings, computing the shifted Laplacian with a multigrid cycle.');
else
    tableCaption = strcat('Model problem 3 with GMRES(',num2str(restart),'): Number of GMRES iterations, matrix-vector products and timings, computing the shifted Laplacian with a multigrid cycle.');
end
tableLabel   = strcat('mg_','restart',restt,'_',dimn,'D_',bdc);

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table1  = {'\begin{table}[t]';...
    '\centering';header;'\toprule'};

row1   = {'$k$ & \multicolumn{3}{c}{CSL(MG)} & \multicolumn{3}{c}{FP(1)+CSL(MG)} & ',...
      '\multicolumn{3}{c}{FP(2)+CSL(MG)}\\'};
row2  = {'& Iter. & MV & Time(s) &  Iter. & MV & Time(s) & ', ...
    '       Iter. & MV & Time(s) \\ \midrule'};

table1 = [table1;row1';row2'];

for i=1:numrows
    for j =1:numcols
        dataValue = num2str(table1_data(i,j),dataFormat{j});
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    table1(end+1) = {[rowStr,' \\']};
end

footer = {'\end{tabular}';['\caption{',tableCaption,'}']; ...
    ['\label{table:',tableLabel,'}'];'\end{table}'};
table1 = [table1;'\bottomrule';footer];

% Save the table to a file in folder ../../../tex_files/tables
namefile = strcat('table_mg_','restart',restt,'_',dimn,'D_',bdc,'.tex');
currentpath = mfilename('fullpath');
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
f = fullfile(currentpath,'tex_files','tables',namefile);

if ( save_flag == 1 )
    fid = fopen(f,'w');
    [nrows,ncols] = size(table1);
    for row = 1:nrows
        fprintf(fid,'%s\n',table1{row,:});
    end
    fclose(fid);
end

% Experiments with the CSL-preconditioned 1D Helmholtz equation with Dirichlet.
% We compare timing and the iteration numbers for solving
%    A M^{-1}x = b
% with gmres, for different implementations of A*M^{-1}*x:
%   - lu: A*(U\(L\x))
%   - ilu: A*(iU\(iL\x))
%   - multi-grid
%
% Further, we do the same when using polynomial acceleration with a truncated
% Faber series s_n(z) of 1/z on a bw set. Then the system to solve is
%   AM^{-1}s_n(AM^{-1}) x = b.
%
% OBSERVATIONS: (Dirichlet boundary conditions)
%   acceleration).
% - With truncated Faber series less iterations are needed.
% - With truncated Faber series less time is needed starting at k = 50.

clear all;
close all;
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

%% Setup parameters
%Setup list of wavenumbers
%wavenum = [10:10:90, 100:100:1500];
%wavenum = [20:20:100,100:100:400];

wavenum = [20:20:100,120,150,200,400,800];
%wavenum = 20:20:100;

degree  = 1:3;

%wavenum  = [20:20:100];
%wavenum  = 50; %% run this when testing changes in the code

%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'dir';
ppw = 15;   %number of points per wavelength
warning off

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 250;

numruns = 2; %number of runs for every experiment (to average later)

%memory allocation for time, residuals
num_dofs    = zeros(length(degree),length(wavenum));
time_lu     = zeros(length(wavenum),1);
time_lu_FS  = zeros(length(degree),length(wavenum));
time_mg     = zeros(length(wavenum),1);
time_mg_FS  = zeros(length(degree),length(wavenum));

iter_lu     = zeros(length(wavenum),2);
iter_lu_FS  = zeros(length(degree),length(wavenum),2);
iter_mg     = zeros(length(wavenum),2);
iter_mg_FS  = zeros(length(degree),length(wavenum),2);

for kk = 1:length(wavenum)
    
    %wavenumber
    k = wavenum(kk);
    factoreps = 0.5;
    poweps = 2;
    eps  = factoreps*k^poweps;
    ppw  = 15; %Choose npf approx ceil(k^3/2) (Pollution free grid)
    
    %number of points in finest grid and number of levels
    [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
    
    %1D
    dim  = 1;
    A    = helmholtz(k,0,npf,bc);
    M    = helmholtz(k,eps,npf,bc);
    op_type = 'gal';
    [mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
    
    b  = zeros(length(M),1); ind = floor(length(M)/2);  b(ind)=1;
    x0 = zeros(size(b));
    
    %multigrid parameters
    npre = 1; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 1;
    
    %Setup LU and MG preconditioners
    tic
    [L, U] = lu(M);     % LU = M
    time_setup_lu=toc;
    
    f_lu  = @(x) A*(U\(L\x));
    f_mg  = @(v) A*feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,1);
    
    %% Polynomial acceleration with truncated Faber series  
    deg = 1; % degree of truncated Faber series of 1/z.
    
    %Parameters for the bratwurst shaped set
    lambda    = -1;     % so that 0 is not in the bw set.
    phi       = pi/2;   % or phi = 0.1*pi; ??
    eps_thick = 0.005;
    [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
    
    %Setup LU_FS, iLU_FS and MG_FS preconditioners  
    for t=1:numruns
        % with f_lu
        tic
        [~, ~, ~, iter_lu(kk,:), resvec_lu] = gmres(f_lu,b,restart,tol,maxit);
        time_lu(kk,1) = time_lu(kk,1)+toc;
        
        % with f_mg
        tic
        [~,~,~,iter_mg(kk,:),resvec_mg] = gmres(f_mg,b,restart,tol,maxit);
        time_mg(kk) = time_mg(kk)+ toc;
        
        for s=1:length(degree)
            deg = degree(s);
            
            f_lu_FS  = @(x)truncFS(f_lu,f_lu(x),deg,M_bw,N_bw,'fun');
            f_mg_FS  = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');
            
            % with f_lu_FS
            tic
            [~,~,~,iter_lu_FS(deg,kk,:),resvec_lu_FS] = gmres(f_lu_FS,b,restart,tol,maxit);
            time_lu_FS(deg,kk,1) = time_lu_FS(deg,kk,1)+ toc;
            
            
            % with f_mg_FS
            tic
            profile on
            [X_mg_FS,~,~,iter_mg_FS(deg,kk,:),resvec_mg_FS] = gmres(f_mg_FS,b,restart,tol,maxit);
            time_mg_FS(deg,kk,1) = time_mg_FS(deg,kk,1)+ toc;
            profile off
        end
        
    end
    
end % of going through different wavenumbers

time_lu     = time_lu/numruns + time_setup_lu;
time_mg     = time_mg/numruns;
time_lu_FS  = time_lu_FS/numruns + time_setup_lu;
time_mg_FS  = time_mg_FS/numruns;

%time_ilu    = time_ilu/numruns + time_setup_ilu;
%time_ilu_FS = time_ilu_FS/numruns + time_setup_ilu;


%% Count matrix-vector operations
% For GMRES on A*x = b, each step requires one multiplication by the matrix A.
% For GMRES on p(A)A*x = b, each step requires one multiplication by the matrix
% p(A)A, i.e., 1 + deg(p) multiplications by A.
% Applying this to the matrix M^{-1} A in our problem, we see that LU, iLU, MG
% have ITER multiplications by M^{-1}A, i.e., M^{-1} is applied ITER times.
% For LU+FS, iLU+FS, MG+FS, we have ITER * (1+d) multiplications by M^{-1}A and
% applications of M^{-1}.

if (isempty(restart))
    rest = 0;
else
    rest = restart;
end

totiter_mg = (iter_mg(:,1)-1)*rest + iter_mg(:,2);
mvop_mg   = totiter_mg;

totiter_lu = (iter_lu(:,1)-1)*rest + iter_lu(:,2);
mvop_lu   = totiter_lu;

totiter_mg_FS = (iter_mg_FS(:,:,1)-1)*rest + iter_mg_FS(:,:,2);
mvop_mg_FS    = diag((degree+1))*totiter_mg_FS;

totiter_lu_FS = (iter_lu_FS(:,:,1)-1)*rest + iter_lu_FS(:,:,2);
mvop_lu_FS    = diag((degree+1))*totiter_lu_FS;

%mvop_lu     = iter_lu(:,2);
%mvop_lu_FS  = (1 + deg)*iter_lu_FS(:,2); %fix this for restarted gmres
%mvop_mg_FS  = (1 + deg)*iter_mg_FS(:,2); %fix this for restarted gmres


%% Plot timings and save the plot in tikz (.tex) and .eps formats
% figure(1)
% plot(wavenum, time_lu, 'k--','linewidth',2)
% plot(wavenum, time_mg, 'k-','linewidth',2)
% plot(wavenum, time_lu_FS(1,:), 'b--','linewidth',2)
% plot(wavenum, time_mg_FS(1,:), 'b-','linewidth',2)
% plot(wavenum, time_lu_FS(2,:), 'r--','linewidth',2)
% plot(wavenum, time_mg_FS(2,:), 'r-','linewidth',2)
% plot(wavenum, time_lu_FS(3,:), 'm--','linewidth',2)
% plot(wavenum, time_mg_FS(3,:), 'm-','linewidth',2)
%
% hold off
% ylabel('Time (s)')
% xlabel('Wavenumber')
% legend('CSL (LU)',...
%        'CSL (MG)', ...
%        'FS + CSL (LU), deg=1', ...
%        'FS + CSL (MG), deg=1',...
%        'FS + CSL (LU), deg=2', ...
%        'FS + CSL (MG), deg=2',...
%        'FS + CSL (LU), deg=3', ...
%        'FS + CSL (MG), deg=3',...
%        'Location','NorthWest')


figure(1)
plot(wavenum, time_lu, 'k-','linewidth',2)
hold on
plot(wavenum, time_lu_FS(1,:), 'b-','linewidth',2)
plot(wavenum, time_lu_FS(2,:), 'r-','linewidth',2)
plot(wavenum, time_lu_FS(3,:), 'm-','linewidth',2)
hold off

ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL (LU)', ...
    'FP(1) + CSL (LU)',...
    'FP(2) + CSL (LU)',...
    'FP(3) + CSL (LU)',...
    'Location','NorthWest')

FS = 16; %font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)
box on


%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the Matlab path
kmin    = num2str(min(wavenum));
kmax    = num2str(max(wavenum));
dimn    = num2str(dim);
degmin  = num2str(min(degree));
degmax  = num2str(max(degree));
restt   = num2str(rest);
bdc     = num2str(bc);

plot_time_tex = strcat('time_lu_vs_k_',kmin,'to',kmax,'_',dimn,'D_',...
    bdc,'deg',degmin,'to',degmax,'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');

%generate file name, for saving in tex_files/figures/new_exp
path = fullfile(currentpath,'..','..','..','..','tex_files','figures','new_exp');

plot_file_tex = fullfile(path,plot_time_tex);

if (save_flag == 1)
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    plot_time_eps = strcat('time_lu_vs_k_',kmin,'to',kmax,'_',dimn,'D_',...
        bdc,'deg',degmin,'to',degmax,'restart',restt,'.eps');
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end

figure(2)
plot(wavenum, time_mg, 'k-','linewidth',2)
hold on
plot(wavenum, time_mg_FS(1,:), 'b-','linewidth',2)
plot(wavenum, time_mg_FS(2,:), 'r-','linewidth',2)
plot(wavenum, time_mg_FS(3,:), 'm-','linewidth',2)
ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL (MG)', ...
    'FP(1) + CSL (MG)',...
    'FP(2) + CSL (MG)',...
    'FP(3) + CSL (MG)',...
    'Location','NorthWest')

FS = 16; %font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)
box on

%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the Matlab path
kmin = num2str(min(wavenum));
kmax = num2str(max(wavenum));
dimn = num2str(dim);
degr = num2str(deg);
bdc  = num2str(bc);

plot_time_tex = strcat('time_mg_vs_k_',kmin,'to',kmax,'_',dimn,'D_',...
    bdc,'deg',degmin,'to',degmax,'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');

%generate file name, for saving in tex_files/figures/new_exp
path = fullfile(currentpath,'..','..','..','..','tex_files','figures','new_exp');

plot_file_tex = fullfile(path,plot_time_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
     plot_time_eps = strcat('time_mg_vs_k_',kmin,'to',kmax,'_',dimn,'D_',...
        bdc,'deg',degmin,'to',degmax,'restart',restt,'.eps');
      
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end


%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(3)
plot(wavenum, totiter_lu(:), 'k-','linewidth',2)
hold on
plot(wavenum, totiter_lu_FS(1,:), 'b-','linewidth',2)
plot(wavenum, totiter_lu_FS(2,:), 'r-','linewidth',2)
plot(wavenum, totiter_lu_FS(3,:), 'm-','linewidth',2)
hold off
ylabel('Number of GMRES iterations')
xlabel('Wavenumber')
legend('CSL (LU)', ...
    'FP(1) + CSL (LU)',...
    'FP(2) + CSL (LU)',...
    'FP(3) + CSL (LU)',...
    'Location','NorthWest')


%title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%num2str(deg)])
FS = 16; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS);
box on

plot_iter_tex = strcat('iter_lu_vs_k_',kmin,'to',kmax,...
    '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
    'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');
%generate file name, for saving in tex_files/figures/new_exp
path          = fullfile(currentpath,'..','..','..','..','tex_files',...
    'figures','new_exp');
plot_file_tex = fullfile(path,plot_iter_tex);

if (save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    plot_iter_eps = strcat('iter_lu_vs_k_',kmin,'to',kmax,...
                      '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
                      'restart',restt,'.eps');
    plot_file_eps = fullfile(path,plot_iter_eps);
    print('-depsc2',plot_file_eps);
end


figure(4)
plot(wavenum, totiter_mg(:), 'k-','linewidth',2)
hold on
plot(wavenum, totiter_mg_FS(1,:), 'b-','linewidth',2)
plot(wavenum, totiter_mg_FS(2,:), 'r-','linewidth',2)
plot(wavenum, totiter_mg_FS(3,:), 'm-','linewidth',2)

hold off
ylabel('Number of GMRES iterations')
xlabel('Wavenumber')
legend('CSL (MG)', ...
    'FP(1) + CSL (MG)',...
    'FP(2) + CSL (MG)',...
    'FP(3) + CSL (MG)',...
    'Location','NorthWest')

%title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
%num2str(deg)])
FS = 16; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS);
box on

plot_iter_tex = strcat('iter_mg_vs_k_',kmin,'to',kmax,...
    '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
    'restart',restt,'.tex');

%get path of current .m file
currentpath  = mfilename('fullpath');
%generate file name, for saving in tex_files/figures/new_exp
path          = fullfile(currentpath,'..','..','..','..','tex_files',...
    'figures','new_exp');
plot_file_tex = fullfile(path,plot_iter_tex);

if (save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', ...
        'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    plot_iter_eps = strcat('iter_mg_vs_k_',kmin,'to',kmax,...
                           '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
                            'restart',restt,'.eps');
   
    plot_file_eps = fullfile(path,plot_iter_eps);
    print('-depsc2',plot_file_eps);
end


% %% Plot relative GMRES residuals vs number of applications of M^{-1}
% % for fixed wavenumber (the last one). For a different wavenumber: grab the
% % vector with GMRES residuals and modify the following plot.
%
% figure(3)
% semilogy(0:mvop_lu, resvec_lu/resvec_lu(1), 'k--')
% hold on
% %semilogy(0:mvop_ilu, resvec_ilu/resvec_ilu(1), 'k-.')
% semilogy(0:mvop_mg, resvec_mg/resvec_mg(1), 'k-')
% semilogy(0:(1+deg):mvop_lu_FS, resvec_lu_FS/resvec_lu_FS(1), 'b--')
% %semilogy(0:(1+deg):mvop_ilu_FS, resvec_ilu_FS/resvec_ilu_FS(1), 'b-.')
% semilogy(0:(1+deg):mvop_mg_FS, resvec_mg_FS/resvec_mg_FS(1), 'b-')
% hold off
% ylabel('Relative residual in GMRES')
% xlabel('Number of applications of M^{-1}')
% legend('Exact inversion of preconditioner', ...
%        'MG preconditioner', ...
%     'Faber series with exact preconditioner', ...
%     'Faber series with MG preconditioner', ...
%     'Location','NorthWest')
% %title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
% %num2str(deg)])
% FS = 14; % font size
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'FontSize',FS)


%% Export .tex table with the timings for several preconditioning methods

%Table with comparison data of LU
table1_data = [wavenum.',totiter_lu,mvop_lu,time_lu,...
               totiter_lu_FS(1,:)',mvop_lu_FS(1,:)',time_lu_FS(1,:)',...
               totiter_lu_FS(2,:)',mvop_lu_FS(2,:)',time_lu_FS(2,:)',...
               totiter_lu_FS(3,:)',mvop_lu_FS(3,:)',time_lu_FS(3,:)'];

numrows      = size(table1_data,1);
numcols      = size(table1_data,2);
tableCaption = strcat('Model problem 1: Number of GMRES iterations, matrix-vector products and timings, computing the shifted Laplacian with an LU factorization.');
tableLabel   = strcat('lu_',dimn,'D_',bdc);

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table1  = {'\begin{table}[ht]';'\resizebox{\textwidth}{!}{';...
    '\centering';header;'\toprule'};

% row1   = {'$k$ & \multicolumn{2}{c}{CSL+LU} & \multicolumn{2}{c}{FS+CSL+LU} & ',...
%           '\multicolumn{2}{c}{CSL+iLU} & \multicolumn{2}{c}{FS+CSL+iLU} & ',...
%     '\multicolumn{2}{c}{CSL+MG} & \multicolumn{2}{c}{FS+CSL+MG}\\'};
% row2  = {'& Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) \\ \midrule'};

row1   = {'$k$ & \multicolumn{3}{c}{CSL(LU)} & \multicolumn{3}{c}{FP(1)+CSL(LU)} & ',...
      '\multicolumn{3}{c}{FP(2)+CSL(LU)}& \multicolumn{3}{c}{FP(3)+CSL(LU)}\\'};
row2  = {'& Iter. & MV & Time(s) &  Iter. & MV & Time(s) & ', ...
    '       Iter. & MV & Time(s) &  Iter. & MV & Time(s) \\ \midrule'};

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

footer = {'\end{tabular}}';['\caption{',tableCaption,'}']; ...
    ['\label{table:',tableLabel,'}'];'\end{table}'};
table1 = [table1;'\bottomrule';footer];

% Save the table to a file in folder ../../../tex_files/tables
namefile = strcat('table_lu_',dimn,'D_',bdc,'.tex');
currentpath = mfilename('fullpath');
f = fullfile(currentpath,'..','..','..','..','tex_files','tables',namefile);

if ( save_flag == 1 )
    fid = fopen(f,'w');
    [nrows,ncols] = size(table1);
    for row = 1:nrows
        fprintf(fid,'%s\n',table1{row,:});
    end
    fclose(fid);
end


%Table with comparison data of MG
table2_data = [wavenum.',totiter_mg,mvop_mg,time_mg,...
               totiter_mg_FS(1,:)',mvop_mg_FS(1,:)',time_mg_FS(1,:)',...
               totiter_mg_FS(2,:)',mvop_mg_FS(2,:)',time_mg_FS(2,:)',...
               totiter_mg_FS(3,:)',mvop_mg_FS(3,:)',time_mg_FS(3,:)'];

numrows = size(table2_data,1);
numcols = size(table2_data,2);
tableCaption = strcat('Model problem 1: Number of GMRES iterations, matrix-vector products and timings, computing the shifted Laplacian with a MG cycle.');
tableLabel   = strcat('mg',dimn,'D_',bdc);

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table2  = {'\begin{table}[t]';'\resizebox{\textwidth}{!}{';...
    '\centering';header;'\toprule'};

row1   = {'$k$ & \multicolumn{3}{c}{CSL(MG)} & \multicolumn{3}{c}{FP(1)+CSL(MG)} & ',...
    '\multicolumn{3}{c}{FP(2)+CSL(MG)} & \multicolumn{3}{c}{FP(3)+CSL(MG)}\\'};
row2   = {' & Iter. & MV & Time(s) & Iter. & MV & Time(s) & ', ...
    '   Iter. & MV & Time(s) & Iter. & MV & Time(s) \\ \midrule'};

table2 = [table2;row1';row2'];

for i=1:numrows
    for j =1:numcols
        dataValue = num2str(table2_data(i,j),dataFormat{j});
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    table2(end+1) = {[rowStr,' \\']};
end

footer = {'\end{tabular}}';['\caption{',tableCaption,'}']; ...
    ['\label{table:',tableLabel,'}'];'\end{table}'};
table2 = [table2;'\bottomrule';footer];

% Save the table to a file in folder ../../../tex_files/tables
namefile = strcat('table_mg_',dimn,'D_',bdc,'.tex');
f = fullfile(currentpath,'..','..','..','..','tex_files','tables',namefile);

if (save_flag == 1)
    fid = fopen(f,'w');
    [nrows,ncols] = size(table2);
    for row = 1:nrows
        fprintf(fid,'%s\n',table2{row,:});
    end
    fclose(fid);
end


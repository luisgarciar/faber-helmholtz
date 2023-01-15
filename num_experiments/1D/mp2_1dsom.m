% CSL-preconditioned 1D Helmholtz with Sommerfeld boundary conditions
%
% We compare timing and the iteration numbers for solving
%    A M^{-1}x = b
% with gmres (M = CSL), for different implementations of A*M^{-1}*x:
%   - lu: A*(U\(L\x))
%   - ilu: A*(iU\(iL\x))
%   - multi-grid
%
% Further, we do the same when using polynomial acceleration with a truncated
% Faber series s_n(z) of 1/z on a bw set. Then the system to solve is
%   AM^{-1}s_n(AM^{-1}) x = b.
%
% OBSERVATIONS: (Sommerfeld boundary conditions)

% For Dirichlet it was:
% - With truncated Faber series less iterations are needed.
% - With truncated Faber series less time is needed starting at k = 50.

save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.
LW = 'LineWidth'; lw = 1;

%% Setup parameters
%Setup list of wavenumbers
wavenum = [20, 40, 60]; % warm up
wavenum = [20:20:100,120,150,200,400,600,800];

degree  = 1:3;

%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'som';
ppw = 15;   %number of points per wavelength

warning off

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 350;

numruns = 1; %number of runs for every experiment (to average later)

%memory allocation for time, residuals
time_lu     = zeros(length(wavenum), 1);
time_lu_FS  = zeros(length(wavenum), length(degree));
time_mg     = zeros(length(wavenum), 1);
time_mg_FS  = zeros(length(wavenum), length(degree));

iter_lu     = zeros(length(wavenum), 2);
iter_lu_FS  = zeros(length(wavenum), length(degree), 2);
iter_mg     = zeros(length(wavenum), 2);
iter_mg_FS  = zeros(length(wavenum), length(degree), 2);

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
            [~,~,~,iter_lu_FS(kk,deg,:),resvec_lu_FS] = gmres(f_lu_FS,b,restart,tol,maxit);
            time_lu_FS(kk,deg,1) = time_lu_FS(kk,deg,1)+ toc;
            
            
            % with f_mg_FS
            tic
            profile on
            [X_mg_FS,~,~,iter_mg_FS(kk,deg,:),resvec_mg_FS] = gmres(f_mg_FS,b,restart,tol,maxit);
            time_mg_FS(kk,deg,1) = time_mg_FS(kk,deg,1)+ toc;
            profile off
        end
        
    end
    
end % of going through different wavenumbers

time_lu     = time_lu/numruns + time_setup_lu;
time_mg     = time_mg/numruns;
time_lu_FS  = time_lu_FS/numruns + time_setup_lu;
time_mg_FS  = time_mg_FS/numruns;


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
mvop_mg_FS    = totiter_mg_FS * diag(degree+1);

totiter_lu_FS = (iter_lu_FS(:,:,1)-1)*rest + iter_lu_FS(:,:,2);
mvop_lu_FS    = totiter_lu_FS * diag(degree+1);


%% Plot timings and save the plot in tikz (.tex) and .eps formats

figure(1)
plot(wavenum, time_lu, 'k-o', LW, lw)
hold on
plot(wavenum, time_lu_FS(:,1), 'b--x', LW, lw)
plot(wavenum, time_lu_FS(:,2), 'r-.+', LW, lw)
plot(wavenum, time_lu_FS(:,3), 'm:s', LW, lw)
hold off

ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL(LU)', ...
    'FP(1)+CSL(LU)',...
    'FP(2)+CSL(LU)',...
    'FP(3)+CSL(LU)',...
    'Location','NorthWest','FontSize',16,'interpreter','latex')

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
% path = fullfile(currentpath,'..','..','..','..','tex_files','figures','new_exp');
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
path = fullfile(currentpath, 'tex_files', 'figures', 'new_exp');

plot_file_tex = fullfile(path,plot_time_tex);

if (save_flag == 1)
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex, 'width', '\figWidth', ...
        'extraaxisoptions',...
        ['xlabel style={font=\scriptsize, at={(axis description cs:0.5,-0.05)},anchor=north},', ...
        'ylabel style={font=\scriptsize, at={(axis description cs:-0.07,.5)},anchor=south},', ...
        'legend style={font=\scriptsize},']);
    
    %Save the plot as .eps
    plot_time_eps = strcat('time_lu_vs_k_',kmin,'to',kmax,'_',dimn,'D_',...
        bdc,'deg',degmin,'to',degmax,'restart',restt,'.eps');
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end

figure(2)
plot(wavenum, time_mg, 'k-o', LW, lw)
hold on
plot(wavenum, time_mg_FS(:,1), 'b--x', LW, lw)
plot(wavenum, time_mg_FS(:,2), 'r-.+', LW, lw)
plot(wavenum, time_mg_FS(:,3), 'm:s', LW, lw)
hold off
ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL(MG)', ...
    'FP(1)+CSL(MG)',...
    'FP(2)+CSL(MG)',...
    'FP(3)+CSL(MG)',...
    'Location','NorthWest','FontSize',16,'interpreter','latex')


%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the matlab path
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
% path = fullfile(currentpath,'..','..','..','..','tex_files','figures','new_exp');
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
path = fullfile(currentpath, 'tex_files', 'figures', 'new_exp');


plot_file_tex = fullfile(path,plot_time_tex);

if ( save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'width', '\figWidth', ...
        'extraaxisoptions',...
        ['xlabel style={font=\scriptsize, at={(axis description cs:0.5,-0.05)},anchor=north},', ...
        'ylabel style={font=\scriptsize, at={(axis description cs:-0.07,.5)},anchor=south},', ...
        'legend style={font=\scriptsize},']);
    
    %Save the plot as .eps
    plot_time_eps = strcat('time_mg_vs_k_',kmin,'to',kmax,'_',dimn,'D_',...
        bdc,'deg',degmin,'to',degmax,'restart',restt,'.eps');
      
    plot_file_eps = fullfile(path,plot_time_eps);
    print('-depsc2',plot_file_eps)
end


%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(3)
plot(wavenum, totiter_lu(:), 'k-o', LW, lw)
hold on
plot(wavenum, totiter_lu_FS(:,1), 'b--x', LW, lw)
plot(wavenum, totiter_lu_FS(:,2), 'r-.+', LW, lw)
plot(wavenum, totiter_lu_FS(:,3), 'm:s', LW, lw)
hold off
ylabel('Number of GMRES iterations')
xlabel('Wavenumber')
legend('CSL(LU)', ...
    'FP(1)+CSL(LU)',...
    'FP(2)+CSL(LU)',...
    'FP(3)+CSL(LU)',...
    'Location','NorthWest','FontSize',16,'interpreter','latex')
ylim([0, 330])

plot_iter_tex = strcat('iter_lu_vs_k_',kmin,'to',kmax,...
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

if (save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex, 'width', '\figWidth', ...
        'extraaxisoptions',...
        ['xlabel style={font=\scriptsize, at={(axis description cs:0.5,-0.05)},anchor=north},', ...
        'ylabel style={font=\scriptsize, at={(axis description cs:-0.09,.5)},anchor=south},', ...
        'legend style={font=\scriptsize},']);
    
    %Save the plot as .eps
    plot_iter_eps = strcat('iter_lu_vs_k_',kmin,'to',kmax,...
                      '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
                      'restart',restt,'.eps');
    plot_file_eps = fullfile(path,plot_iter_eps);
    print('-depsc2',plot_file_eps);
end


figure(4)
plot(wavenum, totiter_mg(:), 'k-o', LW, lw)
hold on
plot(wavenum, totiter_mg_FS(:,1), 'b--x', LW, lw)
plot(wavenum, totiter_mg_FS(:,2), 'r-.+', LW, lw)
plot(wavenum, totiter_mg_FS(:,3), 'm:s', LW, lw)
hold off
ylabel('Number of GMRES iterations')
xlabel('Wavenumber')
legend('CSL(MG)', ...
    'FP(1)+CSL(MG)',...
    'FP(2)+CSL(MG)',...
    'FP(3)+CSL(MG)',...
    'Location','NorthWest','FontSize',16,'interpreter','latex')
ylim([0, 330])

plot_iter_tex = strcat('iter_mg_vs_k_',kmin,'to',kmax,...
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

if (save_flag == 1 )
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex, 'width', '\figWidth', ...
        'extraaxisoptions',...
        ['xlabel style={font=\scriptsize, at={(axis description cs:0.5,-0.05)},anchor=north},', ...
        'ylabel style={font=\scriptsize, at={(axis description cs:-0.09,.5)},anchor=south},', ...
        'legend style={font=\scriptsize},']);
    
    %Save the plot as .eps
    plot_iter_eps = strcat('iter_mg_vs_k_',kmin,'to',kmax,...
                           '_',dimn,'D_',bdc,'deg',degmin,'to',degmax,...
                            'restart',restt,'.eps');
   
    plot_file_eps = fullfile(path,plot_iter_eps);
    print('-depsc2',plot_file_eps);
end


%% Export .tex table with the timings for several preconditioning methods

%Table with comparison data of LU
table1_data = [wavenum.',totiter_lu,mvop_lu,time_lu,...
               totiter_lu_FS(:,1),mvop_lu_FS(:,1),time_lu_FS(:,1),...
               totiter_lu_FS(:,2),mvop_lu_FS(:,2),time_lu_FS(:,2),...
               totiter_lu_FS(:,3),mvop_lu_FS(:,3),time_lu_FS(:,3)];

numrows      = size(table1_data,1);
numcols      = size(table1_data,2);
tableCaption = strcat('Model problem 2: Number of GMRES iterations, matrix-vector products and timings, computing the shifted Laplacian with an LU factorization.');
tableLabel   = strcat('lu_',dimn,'D_',bdc);

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table1  = {'\begin{table}[t]';'\resizebox{\textwidth}{!}{';...
    '\centering';header;'\toprule'};

% row1   = {'$k$ & \multicolumn{2}{c}{CSL+LU} & \multicolumn{2}{c}{FS+CSL+LU} & ',...
%           '\multicolumn{2}{c}{CSL+iLU} & \multicolumn{2}{c}{FS+CSL+iLU} & ',...
%     '\multicolumn{2}{c}{CSL+MG} & \multicolumn{2}{c}{FS+CSL+MG}\\'};
% row2  = {'& Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) \\ \midrule'};

row1   = {'$k$ & \multicolumn{3}{c}{CSL(LU)} & \multicolumn{3}{c}{FP(1)+CSL(LU)} & ',...
      '\multicolumn{3}{c}{FP(2)+CSL(LU)}& \multicolumn{3}{c}{FP(3)+CSL(LU)}\\'};
row2  = {'& Iter. & MV & Time(s) & Iter. & MV & Time(s) & ', ...
    '       Iter. & MV & Time(s) & Iter. & MV & Time(s) \\ \midrule'};

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
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
% path = fullfile(currentpath, 'tex_files', 'figures', 'new_exp');

f = fullfile(currentpath,'tex_files','tables',namefile);

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
               totiter_mg_FS(:,1),mvop_mg_FS(:,1),time_mg_FS(:,1),...
               totiter_mg_FS(:,2),mvop_mg_FS(:,2),time_mg_FS(:,2),...
               totiter_mg_FS(:,3),mvop_mg_FS(:,3),time_mg_FS(:,3)];

numrows = size(table2_data,1);
numcols = size(table2_data,2);
tableCaption = strcat('Model problem 2: Number of GMRES iterations, matrix-vector products and timings, computing the shifted Laplacian with a multigrid cycle.');
tableLabel   = strcat('mg',dimn,'D_',bdc);

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF,...
              '%.0f','%.0f',TF};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table2  = {'\begin{table}[ht]';'\resizebox{\textwidth}{!}{';...
    '\centering';header;'\toprule'};

row1   = {'$k$ & \multicolumn{3}{c}{CSL(MG)} & \multicolumn{3}{c}{FP(1)+CSL(MG)} & ',...
    '\multicolumn{3}{c}{FP(2)+CSL(MG)} & \multicolumn{3}{c}{FP(3)+CSL(MG)}\\'};
row2   = {' & Iter. & MV & Time(s) & Iter. & MV & Time(s) & ', ...
    '         Iter. & MV & Time(s) & Iter. & MV & Time(s) \\ \midrule'};

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

currentpath = mfilename('fullpath');
idcs   = strfind(currentpath,filesep);
currentpath = currentpath(1:idcs(end-3)-1);
f = fullfile(currentpath,'tex_files','tables',namefile);

if (save_flag == 1)
    fid = fopen(f,'w');
    [nrows,ncols] = size(table2);
    for row = 1:nrows
        fprintf(fid,'%s\n',table2{row,:});
    end
    fclose(fid);
end


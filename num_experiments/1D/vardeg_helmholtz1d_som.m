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
% We try different degrees.
%
% OBSERVATIONS:
% - Iterations: Iteration number decreases when applying s_n.
%   Higher n means stronger decrease.
% - Time: MG: Applying s_n is slower than without s_n for k <= 400, and faster
%   for k >= 500.  Larger degrees need more time (but not much).
%   LU and iLU: Needs less time when applying s_n than without, and computation
%   time decreases with incresing deg = 1:3.  (Suprising.)


save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.

%% Setup parameters
%Setup list of wavenumbers
% wavenum = [10:10:90, 100:100:500];
wavenum = [20:20:100, 200:100:1000];
%wavenum = 50; %% run this when testing changes in the code

maxdeg = 2; % Partial sums of degrees 1:maxdeg

%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'som';
ppw = 15;   %number of points per wavelength

warning off

%Parameters for GMRES
restart = [];
tol     = 1e-8;
maxit   = 300;

numruns = 3; %number of runs for every experiment (to average later)

%memory allocation for time, residuals
time_lu     = zeros(length(wavenum),1);
time_lu_FS  = zeros(length(wavenum),maxdeg);
time_ilu    = zeros(length(wavenum),1);
time_ilu_FS = zeros(length(wavenum),maxdeg);
time_mg     = zeros(length(wavenum),1);
time_mg_FS  = zeros(length(wavenum),maxdeg);

iter_lu     = zeros(length(wavenum),2);
iter_lu_FS  = zeros(length(wavenum),maxdeg);
iter_ilu    = zeros(length(wavenum),2);
iter_ilu_FS = zeros(length(wavenum),maxdeg);
iter_mg     = zeros(length(wavenum),2);
iter_mg_FS  = zeros(length(wavenum),maxdeg);

resvec_lu     = zeros(maxit,1);
resvec_ilu    = zeros(maxit,1);
resvec_mg     = zeros(maxit,1);
resvec_lu_FS  = zeros(maxit,1);
resvec_ilu_FS = zeros(maxit,1);
resvec_mg_FS  = zeros(maxit,1);

for kk = 1:length(wavenum)
    %wavenumber
    k = wavenum(kk); eps = 0.5*k^2;
    
    %number of points in finest grid and number of levels
    [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
    
    %1D
    dim  = 1;
    A = helmholtz(k,0,npf,bc);
    M = helmholtz(k,eps,npf,bc);
    op_type = 'gal';
    [mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
    
    b  = zeros(length(M),1); ind = floor(length(M)/2);  b(ind)=1;
    x0 = zeros(size(b));
    
    %multigrid and gmres parameters
    npre = 2; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 1;
    
    %Setup LU, iLU and MG preconditioners
    tic
    [L, U] = lu(M);   % LU = M
    time_setup_lu=toc;
    
    tic
    [iL, iU] = ilu(M);  % iL * iU \approx M.
    time_setup_ilu=toc;
    
    f_lu  = @(x) A*(U\(L\x));
    f_ilu = @(x) A*(iU\(iL\x));
    f_mg  = @(v) A*feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,1);
    
    %% Polynomial acceleration with truncated Faber series
    
    %Parameters for the bratwurst shaped set
    lambda    = -1;     % so that 0 is not in the bw set.
    phi       = pi/2;   % or phi = 0.1*pi; ??
    eps_thick = 0.005;
    [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
    
    for t=1:numruns
        % with f_lu
        tic
        [~, ~, ~, iter_lu(kk,:), resvec_lu] = gmres(f_lu, b, restart, tol, maxit);
        time_lu(kk,1) = time_lu(kk,1)+toc;
        
        % with f_ilu
        tic
        [~, ~, ~, iter_ilu(kk,:), resvec_ilu] = gmres(f_ilu, b, restart, tol, maxit);
        time_ilu(kk,1) = time_ilu(kk,1)+toc;
        
        % with f_mg
        tic
        [~,~,~,iter_mg(kk,:),resvec_mg] = gmres(f_mg, b, restart, tol, maxit);
        time_mg(kk,1) = time_mg(kk,1)+ toc;
        
        for deg = 1:maxdeg; % degree of truncated Faber series of 1/z.
            %Setup LU_FS, iLU_FS and MG_FS preconditioners
            f_lu_FS  = @(x)truncFS(f_lu,f_lu(x),deg,M_bw,N_bw,'fun');
            f_ilu_FS = @(x)truncFS(f_ilu,f_ilu(x),deg,M_bw,N_bw,'fun');
            f_mg_FS  = @(x)truncFS(f_mg,f_mg(x),deg,M_bw,N_bw,'fun');
            
            % with f_lu_FS
            tic
            [~,~,~,itertmp,resvec_lu_FS] = gmres(f_lu_FS,b,restart,tol,maxit);
            time_lu_FS(kk,deg) = time_lu_FS(kk,deg)+ toc;
            iter_lu_FS(kk,deg) = itertmp(1,2);
            
            % with f_ilu_FS
            tic
            [~,~,~,itertmp,resvec_ilu_FS] = gmres(f_ilu_FS,b,restart,tol,maxit);
            time_ilu_FS(kk,deg) = time_ilu_FS(kk,deg)+ toc;
            iter_ilu_FS(kk,deg) = itertmp(1,2);
            
            % with f_mg_FS
            tic
            [X_mg_FS,~,~,itertmp,resvec_mg_FS] = gmres(f_mg_FS,b,restart,tol,maxit);
            time_mg_FS(kk,deg) = time_mg_FS(kk,deg)+ toc;
            iter_mg_FS(kk,deg) = itertmp(1,2);
        end % End of going through different degrees.
    end
    
end % of going through different wavenumbers

time_lu     = time_lu/numruns+time_setup_lu;
time_ilu    = time_ilu/numruns+time_setup_ilu;
time_mg     = time_mg/numruns;
time_lu_FS  = time_lu_FS/numruns+time_setup_lu;
time_ilu_FS = time_ilu_FS/numruns+time_setup_ilu;
time_mg_FS  = time_mg_FS/numruns;

%% Count matrix-vector operations
% For GMRES on A*x = b, each step requires one multiplication by the matrix A.
% For GMRES on p(A)A*x = b, each step requires one multiplication by the matrix
% p(A)A, i.e., 1 + deg(p) multiplications by A.
% Applying this to the matrix M^{-1} A in our problem, we see that LU, iLU, MG
% have ITER multiplications by M^{-1}A, i.e., M^{-1} is applied ITER times.
% For LU+FS, iLU+FS, MG+FS, we have ITER * (1+d) multiplications by M^{-1}A and
% applications of M^{-1}.

mvop_lu  = iter_lu(end,2);
mvop_ilu = iter_ilu(end,2);
mvop_mg  = iter_mg(end,2);
for deg = 1:maxdeg
    mvop_lu_FS(:,deg)  = (1 + deg)*iter_lu_FS(:,deg);
    mvop_ilu_FS(:,deg) = (1 + deg)*iter_ilu_FS(:,deg);
    mvop_mg_FS(:,deg)  = (1 + deg)*iter_mg_FS(:,deg);
end


%% Plot timings and save the plot in tikz (.tex) and .eps formats

LW = 'LineWidth';
lw = 1.5;
FS = 'FontSize';
fs = 22;

figure(1)
plot(wavenum, time_lu, 'k--', LW, lw)
hold on
plot(wavenum, time_ilu, 'k-.', LW, lw)
plot(wavenum, time_mg, 'k-', LW, lw)
plot(wavenum, time_lu_FS, 'b--', LW, lw)
plot(wavenum, time_ilu_FS, 'b-.', LW, lw)
plot(wavenum, time_mg_FS, 'b-', LW, lw)
for deg = 1:maxdeg
    text(wavenum(end), time_lu_FS(end,deg), ['deg=',num2str(deg)], FS, fs)
    text(wavenum(end), time_ilu_FS(end,deg), ['deg=',num2str(deg)], FS, fs)
    text(wavenum(end), time_mg_FS(end,deg), ['deg=',num2str(deg)], FS, fs)
end
hold off
ylabel('Time (s)')
xlabel('Wavenumber')
legend('M^{-1} with LU', 'M^{-1} with iLU', 'M^{-1} with MG', ...
    'Location','NorthWest')
set(gca,'FontSize',fs)


if ( save_flag == 1 )
    %Saving the plot in tikz format
    %Uses matlab2tikz, be sure to include it in the matlab path
    kmin = num2str(min(wavenum));
    kmax = num2str(max(wavenum));
    dimn = num2str(dim);
    degr = num2str(deg);
    bdc  = num2str(bc);
    
    % Get path of current .m file:
    currentpath  = mfilename('fullpath');
    % Generate file name, for saving in tex_files/figures/new_exp
    path = fullfile(currentpath,'..','..','..','tex_files','figures','new_exp');
    
    % File name:
    plot_time_tex = strcat('exp1_time_vs_k_',dimn,'D_',bdc,'_vardeg.tex');
    plot_time_eps = strcat('exp1_time_vs_k_',dimn,'D_',bdc,'_vardeg.eps');
    % File name with path:
    plot_file_tex = fullfile(path,plot_time_tex);
    plot_file_eps = fullfile(path,plot_time_eps);
    
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', 'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    print('-depsc2',plot_file_eps);
end


%% Plot iteration numbers and save as tikz (.tex) and .eps
figure(2)
plot(wavenum, iter_lu(:,2), 'k--', LW, lw)
hold on
plot(wavenum, iter_ilu(:,2), 'k-.', LW, lw)
plot(wavenum, iter_mg(:,2),  'k-', LW, lw)
plot(wavenum, iter_lu_FS, 'b--', LW, lw)
plot(wavenum, iter_ilu_FS, 'b-.', LW, lw)
plot(wavenum, iter_mg_FS, 'b-', LW, lw)
for deg = 1:maxdeg
    text(wavenum(end), iter_lu_FS(end,deg), ['deg=',num2str(deg)], FS, fs)
    text(wavenum(end), iter_ilu_FS(end,deg), ['deg=',num2str(deg)], FS, fs)
    text(wavenum(end), iter_mg_FS(end,deg), ['deg=',num2str(deg)], FS, fs)
end
hold off
ylabel('Number of GMRES iterations')
xlabel('Wavenumber')
legend('M^{-1} with LU', 'M^{-1} with iLU', 'M^{-1} with MG', ...
    'Location','NorthWest')
set(gca,'FontSize',fs)


if ( save_flag == 1 )
    % Get path of current .m file:
    currentpath  = mfilename('fullpath');
    % Generate file name, for saving in tex_files/figures/new_exp
    path = fullfile(currentpath,'..','..','..','tex_files','figures','new_exp');
    
    % File name:
    plot_iter_tex = strcat('exp1_iter_vs_k_',dimn,'D_',bdc,'_vardeg.tex');
    plot_iter_eps = strcat('exp1_iter_vs_k_',dimn,'D_',bdc,'_vardeg.eps');
    % File name with full path:
    plot_file_tex = fullfile(path,plot_iter_tex);
    plot_file_eps = fullfile(path,plot_iter_eps);
    
    %Save as tikz figure in .tex file
    matlab2tikz('filename',plot_file_tex,'standalone',true,'extraaxisoptions',...
        ['xlabel style={font=\Large},', 'ylabel style={font=\Large},']);
    
    %Save the plot as .eps
    print('-depsc2',plot_file_eps);
end


return
%% Export .tex table with the timings for several preconditioning methods

table_data = [wavenum.',iter_lu(:,2),time_lu,iter_lu_FS(:,2),time_lu_FS,...
    iter_ilu(:,2), time_ilu, iter_ilu_FS(:,2), time_ilu_FS, ...
    iter_mg(:,2),time_mg, iter_mg_FS(:,2),time_mg_FS];

numrows = size(table_data,1);
numcols = size(table_data,2);
tableCaption = strcat('Comparison of number of GMRES iterations and timings for 1D Sommerfeld problem');
tableLabel   = strcat(dimn,'D_',bdc);

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f',TF,'%.0f',TF,'%.0f',TF,'%.0f',TF,'%.0f',TF,'%.0f',TF};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table  = {'\begin{table}[h]';'\resizebox{\textwidth}{!}{';...
    '\centering';header;'\toprule'};
row1   = {'$k$ & \multicolumn{2}{c}{CSL+LU} & \multicolumn{2}{c}{FS+CSL+LU} & ',...
    '\multicolumn{2}{c}{CSL+iLU} & \multicolumn{2}{c}{FS+CSL+iLU} & ',...
    '\multicolumn{2}{c}{CSL+MG} & \multicolumn{2}{c}{FS+CSL+MG}\\'};
row2  = {'& Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) &  Iter. & Time(s) \\ \midrule'};
table = [table;row1';row2'];

for i=1:numrows
    for j =1:numcols
        dataValue = num2str(table_data(i,j),dataFormat{j});
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    table(end+1) = {[rowStr,' \\']};
end

footer = {'\end{tabular}}';['\caption{',tableCaption,'}']; ...
    ['\label{table:',tableLabel,'}'];'\end{table}'};
table = [table;'\bottomrule';footer];

% Save the table to a file in folder ../../../tex_files/tables
namefile = strcat('table_',dimn,'D_',bdc,'_vardeg.tex');
currentpath = mfilename('fullpath');
f = fullfile(currentpath,'..','..','..','tex_files','tables',namefile);

if ( save_flag == 1 )
    fid = fopen(f,'w');
    [nrows,ncols] = size(table);
    for row = 1:nrows
        fprintf(fid,'%s\n',table{row,:});
    end
    fclose(fid);
end

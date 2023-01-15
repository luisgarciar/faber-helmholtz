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
save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

%% Setup parameters
%Setup list of wavenumbers
%wavenum  = [20:20:120];

%Parameters of shifted Laplacian
k = 80;
factoreps = 0.6;
poweps = 2;
eps    = factoreps*k^poweps;
npc    = 1;    %number of interior points in the coarsest grid in 1D
ppw    = 12; %pollution free grid
bc     = 'som1';

%Parameters for GMRES
restart = 20;
tol     = 1e-8;
maxit   = 100; % maximum number of gmres cycles
numruns = 1; %number of runs for every experiment (to average later)

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
xlabel('Number of matrix-vector multiplications')
leg{1} = 'CSL';


FS = 16; % font size
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS)

degree = 1:3;
linecolor = {'m','r','b'};

%% Run Experiment
for deg = 1:length(degree)
    
    %Setup of Polynomial acceleration with truncated Faber series
    %Parameters for the bratwurst shaped set
    lambda    = -1; % so that 0 is not in the bw set.
    phi       = pi/2;  % or phi = 0.1*pi; ??
    eps_thick = 0.005;
    [~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
    
    %Setup LU_FS and MG_FS preconditioners
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
        % Applying this to the matrix M^{-1} A in our problem, we see that LU, iLU, MG
        % have ITER multiplications by M^{-1}A, i.e., M^{-1} is applied ITER times.
        % For LU+FS, iLU+FS, MG+FS, we have ITER * (1+d) multiplications by M^{-1}A and
        % applications of M^{-1}.
               
        totiter_mg_FS = (iter_mg_FS(1)-1)*rest + iter_mg_FS(2);
        mvops_mg_FS   = (deg+1)*(0:totiter_mg_FS);
        
        semilogy(mvops_mg_FS, resvec_mg_FS, 'color',linecolor{deg},...
            'linestyle','-','Linewidth',2)
        hold on
        
        leg{deg+1} = strcat('FP+CSL (deg=',num2str(deg),')');
        legend(leg,'Location','NorthEast')
    end
    
end % of going through different wavenumbers


%Saving the plot in tikz format
%Uses matlab2tikz, be sure to include it in the matlab path
wn    = num2str(k);
dimn = num2str(dim);
degmin = num2str(min(degree));
degmax = num2str(max(degree));
restt  = num2str(rest); 
bdc  = num2str(bc);

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

%% Count matrix-vector operations
% For GMRES on A*x = b, each step requires one multiplication by the matrix A.
% For GMRES on p(A)A*x = b, each step requires one multiplication by the matrix
% p(A)A, i.e., 1 + deg(p) multiplications by A.
% Applying this to the matrix M^{-1} A in our problem, we see that LU, iLU, MG
% have ITER multiplications by M^{-1}A, i.e., M^{-1} is applied ITER times.
% For LU+FS, iLU+FS, MG+FS, we have ITER * (1+d) multiplications by M^{-1}A and
% applications of M^{-1}.


% %% Plot relative GMRES residuals vs number of applications of M^{-1}
% % for fixed wavenumber (the last one). For a different wavenumber: grab the
% % vector with GMRES residuals and modify the following plot.
% %
% % figure(3)
% % semilogy(0:mvop_lu, resvec_lu/resvec_lu(1), 'k--')
% % hold on
% % % semilogy(0:mvop_ilu, resvec_ilu/resvec_ilu(1), 'k-.')
% % semilogy(0:mvop_mg, resvec_mg/resvec_mg(1), 'k-')
% % semilogy(0:(1+deg):mvop_lu_FS, resvec_lu_FS/resvec_lu_FS(1), 'b--')
% % % semilogy(0:(1+deg):mvop_ilu_FS, resvec_ilu_FS/resvec_ilu_FS(1), 'b-.')
% % semilogy(0:(1+deg):mvop_mg_FS, resvec_mg_FS/resvec_mg_FS(1), 'b-')
% % hold off
% % ylabel('Relative residual in GMRES')
% % xlabel('Number of applications of M^{-1}')
% % legend('Exact inversion of preconditioner', ...
% %        'MG preconditioner', ...
% %     'Faber series with exact preconditioner', ...
% %     'Faber series with MG preconditioner', ...
% %     'Location','NorthWest')
% % %title(['1D Helmholtz with CSL-preconditioner and poly. accel., deg=',...
% % %num2str(deg)])
% % FS = 14; % font size
% % set(gca,'LooseInset',get(gca,'TightInset'))
% % set(gca,'FontSize',FS)
%
%
% %% Export .tex table with the timings for several preconditioning methods
% table_data = [wavenum.',iter_mg(:,2),time_mg, mvop_mg,...
%     iter_mg_FS(:,2),time_mg_FS,mvop_mg_FS];
%
% numrows = size(table_data,1);
% numcols = size(table_data,2);
% tableCaption = strcat('Comparison of number of GMRES iterations, timings and matrix-times-vector operations for 2D Sommerfeld problem');
% tableLabel   = strcat(dimn,'D_',bdc);
%
% %dataFormat: number of decimals for the variables of each column
% %.0f no decimals, .2f 2 decimals, .3e scientific notation
% TF = '%.3f';    % Time Format.
% dataFormat = {'%.0f','%.0f',TF,'%.0f','%.0f',TF,'%.0f',...
%     '%.0f',TF,'%.0f','%.0f',TF,'%.0f'};
%
% header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
% table  = {'\begin{table}[h]';'\resizebox{\textwidth}{!}{';...
%     '\centering';header;'\toprule'};
%
% row1   = {'$k$ & \multicolumn{3}{c}{CSL+MG} & \multicolumn{3}{c}{FP+CSL+MG} \\ '};
% row2  = {' & Iter. & Time(s) & MatVecs  &  Iter. & Time(s) & Matvecs\\ \midrule'};
%
% table = [table;row1';row2'];
%
% for i=1:numrows
%     for j =1:numcols
%         dataValue = num2str(table_data(i,j),dataFormat{j});
%         if j==1
%             rowStr = dataValue;
%         else
%             rowStr = [rowStr,' & ',dataValue];
%         end
%     end
%     table(end+1) = {[rowStr,' \\']};
% end
%
% footer = {'\end{tabular}}';['\caption{',tableCaption,'}']; ...
%     ['\label{table:',tableLabel,'}'];'\end{table}'};
% table = [table;'\bottomrule';footer];
%
% % Save the table to a file in folder ../../../tex_files/tables
% namefile = strcat('table_',dimn,'D_',bdc,'.tex');
% currentpath = mfilename('fullpath');
% f = fullfile(currentpath,'..','..','..','tex_files','tables',namefile);
%
% if (save_flag == 1)
%     fid = fopen(f,'w');
%     [nrows,ncols] = size(table);
%
%     for row = 1:nrows
%         fprintf(fid,'%s\n',table{row,:});
%     end
%     fclose(fid);
% end


%% Old stuff

% % figure(1)
% % print('-depsc2', ['Helm1D_',bc,'_FS_timing_in_k.eps'])
% % figure(2)
% % print('-depsc2', ['Helm1D_',bc,'_FS_iterations_in_k.eps'])
%
% break

% % Plot the relative residual curves FOR ONE WAVENUMBER k
% figure(1)
% semilogy(0:iter_lu(2), resvec_lu/resvec_lu(1), 'k-')
% hold on
% semilogy(0:iter_ilu(2), resvec_ilu/resvec_ilu(1), 'k-')
% semilogy(0:iter_mg(2), resvec_mg/resvec_mg(1), 'k-')
% semilogy(0:iter_ex_FS(2), resvec_ex_FS/resvec_ex_FS(1), 'b-')
% semilogy(0:iter_lu_FS(2), resvec_lu_FS/resvec_lu_FS(1), 'b-')
% semilogy(0:iter_ilu_FS(2), resvec_ilu_FS/resvec_ilu_FS(1), 'b-')
% semilogy(0:iter_mg_FS(2), resvec_mg_FS/resvec_mg_FS(1), 'b-')
% hold off
% ylabel('relative residual')
% xlabel('iteration')
% legend('f_{ex} = M^{-1}A*x (plain gmres)', 'f_{lu}', 'f_{ilu}', ...
%     'f_{ex} and poly acc with FS', 'f_{lu} and poly acc with FS', ...
%     'f_{ilu} and poly acc with FS')
% title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),')']) FS = 22 ; % font size
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'FontSize',FS);
%
% print('-depsc2', 'Helm1D_rr_FS_different_MinvA.eps')
% table_data = [wavenum', ...,
%               time_lu, iter_lu(:,2),...
%               time_lu_FS, iter_lu_FS(:,2),...
%               time_mg, iter_mg(:,2),...
%               time_mg_FS, iter_mg_FS(:,2)];

% name_table = strcat('table_time_vs_k_',dimn,'D_',...
%                  bdc,'_','kmin',kmin,'_',...
%                  'kmax',kmax,'.tex');
% clear input;
% input.data = table_data;
%
% %We want a complete LaTex document (change this later)
% input.makeCompleteLatexDocument = 1;
% input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
%
% %Set column labels (use empty string for no label):
% input.tableColLabels = {'k','Time Exact CSL (s)','Iter. Exact CSL (s)',...
%                         'Time Faber with Exact CSL (s)','Iter. Faber series with Exact CSL',...
%                         'Time MG CSL (s)','Iter. MG-CSL',...
%                         'Time Faber series with MG-CSL(s)','Iter. Faber series with MG CSL'};
%
% %Set row labels (use empty string for no label):
% %input.tableRowLabels = {'','','','','','','','','',''};
%
% % Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
% input.tableColumnAlignment = 'c';
%
% % Switch table borders on/off (borders are enabled by default):
% input.tableBorders = 1;
%
% % Formatting-string to set the precision of the table values:
% % For using different formats in different rows use a cell array like
% % {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% % where myFormatString_ are formatting-strings and numberOfValues_ are the
% % number of table columns or rows that the preceding formatting-string applies.
% % Please make sure the sum of numberOfValues_ matches the number of columns or
% % rows in input.tableData!
% %
% input.dataFormat = {'%.0f',1,'%.2f',1,'%.0f',1,'%.2f',1,'%.0f',1,'%.2f',1,'%.0f',1,'%.4f',1,'%.0f',1,}; % three digits precision for first two columns, one digit for the last
%
% % Switch table borders on/off (borders are enabled by default):
% input.tableBorders = 1;
%
% %Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
%
% %Generate LaTex code
% latex = latexTable(input);
%
% % save LaTex code as file
% fid=fopen(name_table,'w');
% [nrows,ncols] = size(latex);
% for row = 1:nrows
%     fprintf(fid,'%s\n',latex{row,:});
% end
% fclose(fid);
%
% %fprintf('\n... your LaTex code has been saved as ''Table1.tex'' in your working directory\n');


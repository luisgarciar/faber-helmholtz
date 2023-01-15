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

% Computation: mp3_cluster_iter_time
% Postprocessing: this file

load('mp3_results_cluster')

save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.
LW = 'LineWidth'; lw = 1;

%% Plot timings and save the plot in tikz (.tex) and .eps formats
figure(1)
plot(wavenum, time_mg, 'k-o', LW, lw);
hold on
plot(wavenum, time_mg_FS(:,1), 'b-o', LW, lw);
plot(wavenum, time_mg_FS(:,2), 'r-o', LW, lw);
hold off

ylabel('Time (s)')
xlabel('Wavenumber')
legend('CSL(MG)', ...
    'FP(1)+CSL(MG)',...
    'FP(2)+CSL(MG)',...
    'Location','NorthWest');

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
path          = fullfile(currentpath,'..','..','..','..','tex_files',...
    'figures','new_exp');

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
plot(wavenum, totiter_mg_FS(:,1), 'b-o', LW, lw);
plot(wavenum, totiter_mg_FS(:,2), 'r-o', LW, lw);
hold off

ylabel('Number of iterations')
xlabel('Wavenumber')
legend('CSL(MG)', ...
    'FP(1)+CSL(MG)',...
    'FP(2)+CSL(MG)',...
    'Location','NorthWest')

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
path          = fullfile(currentpath,'..','..','..','..','tex_files',...
    'figures','new_exp');

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
f = fullfile(currentpath,'..','..','..','..','tex_files','tables',namefile);

if ( save_flag == 1 )
    fid = fopen(f,'w');
    [nrows,ncols] = size(table1);
    for row = 1:nrows
        fprintf(fid,'%s\n',table1{row,:});
    end
    fclose(fid);
end

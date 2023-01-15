%test for opening files and saving tables

table1 = rand(10,9);

numrows = size(table1,1);
numcols = size(table1,2);
tableCaption = strcat('test table');
tableLabel   = strcat('lalalabel');

%dataFormat: number of decimals for the variables of each column
%.0f no decimals, .2f 2 decimals, .3e scientific notation
TF = '%.3f';    % Time Format.
dataFormat = {'%.0f','%.0f',TF,'%.0f',...
                     '%.0f',TF,'%.0f',...
                     '%.0f',TF,'%.0f'};

header = ['\begin{tabular}','{',repmat(('c'),1,numcols),'}'];
table  = {'\begin{table}[h]';'\resizebox{\textwidth}{!}{';...
    '\centering';header;'\toprule'};

row1   = {'$0$ & \multicolumn{3}{c}{123} & \multicolumn{3}{c}{456} & \multicolumn{3}{c}{789} \\ '};
row2  =  {' 0 & 1 & 2 & 3  & 4 & 5 & 6 & 7 & 8 & 9\\ \midrule'};

table = [table;row1';row2'];

for i=1:numrows
    for j =1:numcols
        dataValue = num2str(table1(i,j),dataFormat{j});
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    table(end+1) = {[rowStr,'\\']};
end

footer = {'\end{tabular}}';['\caption{',tableCaption,'}']; ...
    ['\label{table:',tableLabel,'}'];'\end{table}'};
table = [table;'\bottomrule';footer];


namefile = strcat('test.tex');
currentpath = mfilename('fullpath');
f = fullfile(currentpath,'..','..','..','..','tex_files','tables',namefile);

[fid,message] = fopen(f,'w');
message;

    [nrows,ncols] = size(table1);
    for row = 1:nrows
        fprintf(fid,'%s\n',table{row,:});
    end
    fclose(fid);

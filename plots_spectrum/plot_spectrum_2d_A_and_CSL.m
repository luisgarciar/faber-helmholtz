% Plot the spectrum of A, A*M^{-1} in 2D.
%
% This is a short version of the file plot_spectrum_2d.m.
% This file is intended for generating the figures as tikzpictures (not
% standalong pictures).
%
% Boundary conditions: Dirichlet, Sommerfeld.

save_flag = 1;  % save_flag=1: save plots and table, =0 do not save.

k = 15;         % 20 works too, but is slow.
eps = 0.5*k^2;  % Shift

%number of interior points in coarsest grid in one dim
npc = 1;
bvc = {'dir'};
% bc = 'som';
ppw = 10;       %number of points per wavelength

% Number of points in finest grid and number of levels:
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);

for jj = 1:length(bvc)
    bc = bvc{jj};   % run through both boundary conditions
    
    % Matrices (2D):
    A = helmholtz2(k,0,npf,npf,bc);
    M = helmholtz2(k,eps,npf,npf,bc);
    AMinv = full(M\A);
    
    % Eigenvalues:
    eigA = eig(full(A));
    eigAM = eig(AMinv);
    
    % Plot:
    figure(1)
    plot(real(eigA), imag(eigA), 'b*','MarkerSize',12)
    axis equal
    grid on
    
    figure(2)
    plot(real(eigAM), imag(eigAM), 'b*','MarkerSize',12)
    axis([-0.2 1.2 -0.7 0.7])
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    set(gca,'Xtick',[0 0.5 1], 'Ytick',[-0.5 0. 0.5], 'Fontsize',15);
    
    if (save_flag == 1 )
        % Build path:
        currentpath = mfilename('fullpath'); % path of current m-file
        mypath = fullfile(currentpath,'..','..','tex_files','figures');
        
        %         figure(1)
        %         % Save as eps:
        %         filename = strcat('spectrum_Helm2d_',bc,'.eps');
        %         print('-depsc2', fullfile(mypath,filename))
        %         % Save as tikzpicture (in .tex file):
        %         filename = strcat('spectrum_Helm2d_',bc,'.tex');
        %         matlab2tikz('filename', fullfile(mypath,filename), ...
        %             'width', '\figWidth')
        
        figure(2)
        % Save as eps:
        %         filename = strcat('spectrum_CSL2d_',bc,'.eps');
        %         print('-depsc2', fullfile(mypath,filename))
        % Save as tikzpicture (in .tex file):
        filename = strcat('spectrum_CSL2d_',bc,'.tex');
        matlab2tikz('filename', fullfile(mypath,filename), ...
            'width', '\figWidth')
        % Note: the rectangle command does not translate to a circle in
        % matlab2tikz, so replace manually the corresponding line by
        %
        % \draw[] (0.5, 0) circle (0.5);
    end
    
end








return
%% Plots only eig(A), eig(M^{-1}A)
% figure(1)
% plot(real(eigA), imag(eigA), 'ko')
% grid on
% %axis([0 1 -0.5 0.5])
% xlabel('real(\lambda)','FontSize',20)
% ylabel('imag(\lambda)','FontSize',20)
% set(gca,'FontSize',20)

figure(2)
plot(real(eigAM), imag(eigAM), 'ko')
axis equal
axis([0 1 -0.5 0.5])
set(gca,'Xtick',[0 0.5 1]);
set(gca,'Ytick',[-0.5 0. 0.5]);
xlabel('real(\lambda)','FontSize',20)
ylabel('imag(\lambda)','FontSize',20)
%title('Spectrum of A')
set(gca,'FontSize',20)
rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])



return
%% Old save figures


if ( save_flag == 1 )
    % Build path:
    currentpath = mfilename('fullpath'); % path of current m-file
    mypath = fullfile(currentpath,'..','..','tex_files','figures');
    
    
    figure(1)
    % Save as eps:
    filename = strcat('spectrum_Helm2d_',bc,'.eps');
    print('-depsc2', fullfile(mypath,filename))
    % Save as tikz figure in .tex file
    filename = strcat('spectrum_Helm2d_',bc,'.tex');
    matlab2tikz('filename', fullfile(mypath,filename), 'standalone', true, ...
        'extraaxisoptions',...
        ['xlabel style={font=\Large},', 'ylabel style={font=\Large},']);
    
    % KEEP THIS PART FOR THE ADDITIONAL MATLAB2TIKZ-OPTIONS.
    
    figure(2)
    % Save as eps:
    filename = strcat('spectrum_CSL2d_',bc,'.eps');
    print('-depsc2', fullfile(mypath,filename))
    % Save as tikz figure in .tex file
    filename = strcat('spectrum_CSL2d_',bc,'.tex');
    matlab2tikz('filename', fullfile(mypath,filename), 'standalone', true, ...
        'extraaxisoptions',...
        ['xlabel style={font=\Large},', 'ylabel style={font=\Large},']);
end


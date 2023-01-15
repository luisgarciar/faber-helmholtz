% Plot the spectrum of M\A in 1D.
% 
% Boundary conditions: Dirichlet, Sommerfeld.

wavenum = 80;
bc = 'dir';
% bc = 'som';
save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.

%number of interior points in coarsest grid in one dim
npc = 1;
ppw = 15; %number of points per wavelength

for kk = 1:length(wavenum)
    k = wavenum(kk);    % Wavenumber
    eps = 0.5*k^2;     % Shift
    
    % Number of points in finest grid and number of levels:
    [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
    
    % Matrices (1D):
    A = helmholtz(k,0,npf,bc);
    M = helmholtz(k,eps,npf,bc);
    
    % Eigenvalues:
    eigval = eig(full(M\A));
    
    % Plot:
    figure()
    plot(real(eigval), imag(eigval), 'ko')
    title(['k = ',num2str(k)])
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    largefiglabels
end

if ( save_flag == 1 )
    % Build path:
    currentpath = mfilename('fullpath'); % path of current m-file
    mypath = fullfile(currentpath,'..','..','tex_files','figures');
    
    % Save as eps:
    filename = strcat('spectrum_CSL1d_',bc,'.eps');
    print('-depsc2', fullfile(mypath,filename))
    
    % Save as tikz figure in .tex file
    filename = strcat('spectrum_CSL1d_',bc,'.tex');
    matlab2tikz('filename', fullfile(mypath,filename), 'standalone', true, ...
        'extraaxisoptions',...
        ['xlabel style={font=\Large},', 'ylabel style={font=\Large},']);
end

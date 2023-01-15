% Plot the spectrum of M\A in 2D, and with inclusion set.
%
% Boundary conditions: Dirichlet, Sommerfeld.

wavenum = 15;   % 20 works too, but is slow.
save_flag = 0;  % save_flag=1: save plots and table,=0 do not save.

%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'dir';
ppw = 12;   %number of points per wavelength

%Parameters for the bratwurst shaped set
lambda    = -1; % so that 0 is not in the bw set.
phi       = pi/2;  % or phi = 0.1*pi; ??
eps_thick = 0.02;
[psi, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, eps_thick);
ucirc = chebfun(@(t) exp(1i*t), [0, 2*pi]);
bw = (psi(ucirc) + 1)/2;

for kk = 1:length(wavenum)
    k = wavenum(kk);    % Wavenumber
    eps = 0.5*k^2;      % Shift
    
    % Number of points in finest grid and number of levels:
    [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
    
    % Matrices (2D):
    A = helmholtz2(k,0,npf,npf,bc);
    M = helmholtz2(k,eps,npf,npf,bc);
    
    % Eigenvalues:
%     eigA = eig(full(A));
    eigAM = eig(full(M\A));
        
    figure(1)
    plot(real(eigAM), imag(eigAM), 'k*','MarkerSize',15)
   axis([-0.1 1.1 -0.6 0.6])
   axis equal
    set(gca,'Xtick',[0 0.5 1]);
    set(gca,'Ytick',[-0.5 0. 0.5]);
    xlabel('Real(z)','FontSize',20)
    ylabel('Imag(z)','FontSize',20)
    set(gca,'FontSize',20)
   

    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    hold on
    plot(bw, 'b-','LineWidth',3)
    hold off
    %largefiglabels
end

if ( save_flag == 1 )
    % Build path:
    currentpath = mfilename('fullpath'); % path of current m-file
    mypath = fullfile(currentpath,'..','..','tex_files','figures');
    
    figure(1)
    % Save as eps:
    filename = strcat('spectrum_incluset_CSL2d_',bc,'.eps');
    print('-depsc2', fullfile(mypath,filename))
    % Save as tikz figure in .tex file
    filename = strcat('spectrum_incluset_CSL2d_',bc,'.tex');
    matlab2tikz('filename', fullfile(mypath,filename), 'standalone', true, ...
        'extraaxisoptions',...
        ['xlabel style={font=\Large},', 'ylabel style={font=\Large},']);
end

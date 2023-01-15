% Plot the spectrum of A, A*M^{-1}, and s_n(A M^{-1}) * A M^{-1} in 2D.
%
% Boundary conditions: Dirichlet, Sommerfeld.

wavenum = 40;   % 20 works too, but is slow.
save_flag = 0;  % save_flag=1: save plots and table, =0 do not save.

%number of interior points in coarsest grid in one dim
npc = 1;
bc = 'som';
ppw = 10;   %number of points per wavelength

%Parameters for the bratwurst shaped set
lambda    = -1; % so that 0 is not in the bw set.
phi       = pi/2;  % or phi = 0.1*pi; ??
eps_thick = 0.005;
[psi, ~, ~, Mbw, Nbw] = bw_map(lambda, phi, eps_thick);
% ucirc = chebfun(@(t) exp(1i*t), [0, 2*pi]);
% bw = (psi(ucirc) + 1)/2;

for kk = 1:length(wavenum)
    k = wavenum(kk);    % Wavenumber
    eps = 0.5*k^2;      % Shift
    
    % Number of points in finest grid and number of levels:
    [npf,numlev] = fd_npc_to_npf(npc,k,ppw);
    
    % Matrices (2D):
    A = helmholtz2(k,0,npf,npf,bc);
    M = helmholtz2(k,eps,npf,npf,bc);
    
    AMinv = full(M\A);
    s1AM = truncFS(AMinv, AMinv, 1, Mbw, Nbw, 'mat');   % s_1(AM^{-1})*AM^{-1}
    s2AM = truncFS(AMinv, AMinv, 2, Mbw, Nbw, 'mat');   % s_2(AM^{-1})*AM^{-1}
    s3AM = truncFS(AMinv, AMinv, 3, Mbw, Nbw, 'mat');
    s4AM = truncFS(AMinv, AMinv, 4, Mbw, Nbw, 'mat');
    
    % Eigenvalues:
    eigA = eig(full(A));
    eigAM = eig(AMinv);
    eigs1AM = eig(s1AM);
    eigs2AM = eig(s2AM);
    eigs3AM = eig(s3AM);
    eigs4AM = eig(s4AM);
    
    %%
    figure(1)
    subplot(2,3,1)
    plot(real(eigA), imag(eigA), 'k+')
    title('Spectrum of A')
    grid on
    
    subplot(2,3,2)
    plot(real(eigAM), imag(eigAM), 'k+')
    title('Spectrum of M^{-1}A')
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    axis equal
    
    subplot(2,3,3)
    plot(real(eigs1AM), imag(eigs1AM), 'k+')
    title('Spectrum of s_1(M^{-1}A)M^{-1}A')
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    
    subplot(2,3,4)
    plot(real(eigs2AM), imag(eigs2AM), 'k+')
    title('Spectrum of s_2(M^{-1}A)M^{-1}A')
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    
    subplot(2,3,5)
    plot(real(eigs3AM), imag(eigs3AM), 'k+')
    title('Spectrum of s_3(M^{-1}A)M^{-1}A')
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    
    subplot(2,3,6)
    plot(real(eigs4AM), imag(eigs4AM), 'k+')
    title('Spectrum of s_4(M^{-1}A)M^{-1}A')
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    
%     print('-depsc2', '../tex_files/figures/spectra_helmholtz2d_som_k15.eps')
    
    %return
    % Plot:
    figure(1)
    plot(real(eigA), imag(eigA), 'ko')
    grid on
    largefiglabels
    
    figure(2)
    plot(real(eigAM), imag(eigAM), 'ko')
    axis equal
    rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    largefiglabels
end

%% Plots only eig(A), eig(M^{-1}A)
  figure(1)
    plot(real(eigA), imag(eigA), 'ko')
    grid on
    %axis([0 1 -0.5 0.5])
    xlabel('real(\lambda)','FontSize',20)
    ylabel('imag(\lambda)','FontSize',20)
    set(gca,'FontSize',20)
    
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

    %rectangle('Position', [0, -0.5, 1, 1], 'Curvature', [1 1])
    
    
%%
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


%% Plots for Presentation Cosse
    figure(4)
    subplot(1,4,1)
    plot(real(eigAM), imag(eigAM), 'k+')
    title('Spectrum of M^{-1}A')
    axis([0 1 -0.5 0.5])
   % rectangle('Position', [0, -0.5, -0.5,0.5], 'Curvature', [1 1])
    axis equal
    
    subplot(1,4,2)
    plot(real(eigs1AM), imag(eigs1AM), 'k+')
    title('Spectrum of s_1(M^{-1}A)M^{-1}A')
    axis([0 1 -0.5 0.5])
    axis equal
   % rectangle('Position', [0, -0.5, -0.5, 0.5], 'Curvature', [1 1])
    
    subplot(1,4,3)
    plot(real(eigs2AM), imag(eigs2AM), 'k+')
    title('Spectrum of s_2(M^{-1}A)M^{-1}A')
    axis([0 1 -0.5 0.5])
    axis equal
    %rectangle('Position', [0, -0.5, -0.5, 0.5], 'Curvature', [1 1])
    
    subplot(1,4,4)
    plot(real(eigs3AM), imag(eigs3AM), 'k+')
    title('Spectrum of s_3(M^{-1}A)M^{-1}A')
    axis([0 1 -0.5 0.5])
    axis equal
    %rectangle('Position', [0, -0.5, -0.5, 0.5], 'Curvature', [1 1])
        
    
   

% 
% Experiment:
% Approximate 1/z by a truncated Faber series on a bw set.
% 

save_flag = 0;

% Parameters for the bratwurst shaped set
lambda = -1;
phi = pi/2;
eps_thick = 0.005;
[psi, ~, capacity, M, N] = bw_map(lambda, phi, eps_thick) ;

kmax = 35 ;
F = get_Fpoly(kmax, M, N) ; %% this gets the actual polynomials (in z)

% Control figure  --  visualize the set
FS = 22 ; %% font size

n = 2^12 ;
% unit_circle = exp(1i * (0: 2*pi/n : 2*pi - 2*pi/n).' ) ;  %% column !
unit_circle = exp(2i*pi*(0:n-1)/n).';
circle = unit_circle/2 + 1/2;
bdry_E = (psi(unit_circle) + 1)/2 ;

figure(1)
plot(real(bdry_E), imag(bdry_E), 'k.')
hold on
plot( real(circle), imag(circle), 'k--')
hold off
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS);



% Compute partial sums of Faber series of 1/z

alpha = - N - sqrt(N^2-1) ;
rho = abs(alpha) ;

sn = 0;
for k=1:kmax+1
    sn = polyadd(sn, (M+alpha)/(alpha^k) * F{k}); % a_{k-1} * F_{k-1} 
    err(k) = max(abs( polyval(sn, bdry_E) - 1./bdry_E )) ;
end

figure(2)
ph_fs = semilogy(0:kmax, err, 'k-') ;
hold on
ph_rho = plot(0:kmax, (1/rho).^(0:kmax), 'k--') ;
hold off
h = legend([ph_fs, ph_rho], '$\Vert z^{-1} - s_n(z) \Vert_E$', ...
    '$\rho^{-k}$');
set(h, 'Interpreter', 'Latex')
% set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'FontSize',FS);

if ( save_flag == 1 )
    print('-dpng', 'approx_err.png')
end

%% err < 1e-2:

I = find(err < 1e-2, 1)



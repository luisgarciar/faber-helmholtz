%
% Zeros of truncated Faber series of 1/z
%

% Parameters for the bratwurst shaped set
lambda    = -1;     % so that 0 is not in the inclusion set.
phi       = pi/2;   % or phi = 0.1*pi; ??
eps_thick = 0.005;
[psi, ~, ~, M, N] = bw_map(lambda, phi, eps_thick);
rho = N + sqrt(N^2 - 1);
S = (M*N-1)/(N-M);

% Get partial sums s1 and s2:
s1 = scalar_truncFS(1, M, N);
s2 = scalar_truncFS(2, M, N);

% Get zeros:
zer_s1 = roots(s1)
zer_s2 = roots(s2)

% For s1: compare with analytic expression of zero:
zer_s1_dir = (2*N + rho + S)/(2*(N-M));
err_s1_zer = zer_s1 - zer_s1_dir;


%% Get partial sums s1 and s2 via their coefficients.
% This was for checking accuracy of expressions for s1 and s2.
% s1_dir = -(rho-M)/(rho^2) * [ 2*(N-M), - 2*N - S - rho];
% a = (2*(N-M))^2;
% b = - 8*N^2 + 4*M*N + 4 - 2*rho*(N-M);
% c = 4*N^2 - 2 - S^2 + 2*N*rho + S*rho + rho^2;
% s2_dir = ((rho-M)/(rho^3)) * [a, b, c];
% % Check that both computations give same polynomial:
% err_s1 = norm(s1 - s1_dir, inf)
% err_s2 = norm(s2 - s2_dir, inf)


%% Plotting:

% Get 'bratwurst' contour:
npts = 2^10;
ucirc = exp(2i*pi*(1:npts)/npts);
bw = (psi(ucirc) + 1)/2;

figure()
plot((ucirc+1)/2, 'k--')
hold on
plot(bw, 'k-')
plot(real(zer_s1), imag(zer_s1), 'bo')
plot(real(zer_s2), imag(zer_s2), 'ro')
hold off
axis equal, grid on


%Zeros of truncated series of 1/z
%(sigma = 1)
% We study the dependence of the zeros on the
% parameter N = sec(phi/4) > 1
%
%

reLambda = @(N) (5*N.^2 + N.*sqrt(N.^2 -1)-2)./(4*N.^2);
imLambda = @(N) sqrt((6*N.^4 + 6*N.^3.*sqrt(N.^2 -1) + 5.*N.^2 - 8))./(4*N.^2);

Nstar = 10000; 
N     = linspace(1,Nstar, 1000);

x = reLambda(N); y=imLambda(N);

figure(2)
plot((ucirc+1)/2, 'k--')
hold on
plot(x,y, 'b.')
plot(x,-y, 'r.')
hold off
axis equal, grid on






%







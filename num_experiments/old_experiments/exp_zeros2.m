
% Plot the disk |w-0.5| <= 0.5
npts = 2^10;
ucirc = exp(2i*pi*(1:npts)/npts);

figure(1)
plot((ucirc+1)/2, 'k--')
axis equal, grid on


% Parameters for the bratwurst shaped set
lambda    = -1;     % so that 0 is not in the inclusion set.
phi       = pi/2;   % or phi = 0.1*pi; ??
P = tan(phi/4) + 1/cos(phi/4);
eps_max = tan(phi/4) * (1 + tan(phi/8));    % max thickness

% Different thickness:
epsvals = linspace(0, eps_max, 10);
epsvals(end) = [];      % eps_max is not allowed

for jj = 1:length(epsvals)
    eps_thick = epsvals(jj);
    
    % Compute M, N:
    sigma = 1 + eps_thick ;
    N = (P/sigma + sigma/P)/2;  %% N_eps
    M = (sigma^2 - 1)/(2*sigma * tan(phi/4)); %% M_eps
    
    rho = N + sqrt(N^2 - 1);
    S = (M*N - 1)/(N - M);
    
    % Get partial sums s1 and s2:
    s1 = scalar_truncFS(1, M, N);
    s2 = scalar_truncFS(2, M, N);
    
    a = 4*(N-M)^2;
    b = -8*N^2 + 4*M*N + 4 - 2*rho*(N-M);
    c = 4*N^2 - 2 - S^2 + 2*N*rho + rho*S + rho^2;
    discr(jj,1) = (b^2 - 4*a*c)/4;
    nosqrt(jj,1) = 8 + 11*M^2 - 22*M*N - 5*N^2 - 14*M^2*N^2 + 28*M*N^3 - 6*N^4;
    
%     assert(discr < 0)
    err(jj,1) = 8 + 11*M^2 - 22*M*N - 5*N^2 - 14*M^2*N^2 + 28*M*N^3 - 6*N^4 ...
        - 6*N*(N-M)^2 * sqrt(N^2 - 1) - discr(jj,1);
    
    % Get zeros:
    zer_s1 = roots(s1);
    zer_s2 = roots(s2);
    
    % For s1: compare with analytic expression of zero:
    zer_s1_dir = (2*N + rho + S)/(2*(N-M));
    err_s1_zer = zer_s1 - zer_s1_dir;
    
    figure(1)
    hold on
    plot(real(zer_s1), imag(zer_s1), 'bo')
    plot(real(zer_s2), imag(zer_s2), 'ro')
    hold off
    drawnow
%     pause(0.3)
end

8 + 11*M^2 - 22*M*N - 5*N^2 - 14*M^2*N^2 + 28*M*N^3 - 6*N^4 ...
    - 6*N*(N-M)^2 * sqrt(N^2 - 1);



return

%% Varying M and N independtly
% is not helpful: the expression for the discriminant becomes positive

% [M, N] = meshgrid(linspace(0, 1, 100), linspace(1, 2, 100));
% F = 8 + 11*M.^2 - 22*M.*N - 5*N.^2 - 14*M.^2.*N.^2 + 28*M.*N.^3 - 6*N.^4 ...
%     - 6*N.*(N-M).^2 .* sqrt(N.^2 - 1);
% min(min(F))
% max(max(F))

%% Fix sigma = 1 (arc) and vary phi:

sigma = 1;  % Fix
phiv = linspace(0, 2*pi, 1002);
phiv(1) = []; 
phiv(end) = [];     % phi in ]0, 2*pi[
F = zeros(size(phiv));

for jj = 1:length(phiv)
    phi = phiv(jj);
    P = tan(phi/4) + 1/cos(phi/4);
    
    N = (P/sigma + sigma/P)/2;  %% N_eps
    M = (sigma^2 - 1)/(2*sigma * tan(phi/4)); %% M_eps
    
    rho = N + sqrt(N^2 - 1);
    S = (M*N - 1)/(N - M);
    
    F(jj) = 8 + 11*M^2 - 22*M*N - 5*N^2 - 14*M^2*N^2 + 28*M*N^3 - 6*N^4; % ...
%     - 6*N.*(N-M).^2 .* sqrt(N.^2 - 1);
end
max(F)
plot(phiv, F)
xlabel('\phi')
% With sigma = 1 fixed, F is decreasing in phi in ]0, 2*pi[.

%% Fix angle and vary sigma in [1, sigma_max[.

phi = 0.9 * pi;     % Fix angle
P = tan(phi/4) + 1/cos(phi/4);
eps_max = tan(phi/4) * (1 + tan(phi/8));    % max thickness
sigma_max = eps_max + 1;

% Different thickness:
sigmav = linspace(1, sigma_max, 100);
% sigmav(end) = [];      % sigma_max is not allowed
F = zeros(size(sigmav));

for jj = 1:length(sigmav)
    sigma = sigmav(jj);
    
    N = (P/sigma + sigma/P)/2;  %% N_eps
    M = (sigma^2 - 1)/(2*sigma * tan(phi/4)); %% M_eps
    
    rho = N + sqrt(N^2 - 1);
    S = (M*N - 1)/(N - M);
    
    F(jj) = 8 + 11*M^2 - 22*M*N - 5*N^2 - 14*M^2*N^2 + 28*M*N^3 - 6*N^4; % ...
%     - 6*N.*(N-M).^2 .* sqrt(N.^2 - 1);
end

max(F)
plot(sigmav, F)
xlabel('\sigma')


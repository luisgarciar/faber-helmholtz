
% Parameters of bratwurst set in article:
lambda = -1;
phi = pi/2;
sigma = 1.005;

% get M, N
[~, ~, ~, M_bw, N_bw] = bw_map(lambda, phi, sigma-1);

s1 = scalar_truncFS(1, M_bw, N_bw)
s2 = scalar_truncFS(2, M_bw, N_bw)
s3 = scalar_truncFS(3, M_bw, N_bw)

figure(1)
plotcp(roots(s3), 'x')
shg
rectangle('Position', [0,-0.5,1,1], 'Curvature', [1, 1])
axis equal
grid on

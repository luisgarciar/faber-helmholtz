function [N,iter,time,err,errHist] = amgtest(node,elem,modification,option,varargin)
%% AMGTEST test algebraic mutligrid solver
%
% amgtest(node,elem) test AMG solver for solving Poisson equation $-\Delta u = 1$
% with zero Dirichlet boundary condition.
%
% amgtest(node,elem,modification) test more cases:
% - 0: Neumann problem
% - 1: Dirichlet problem (default choice)
% - 2: time discretization added
% - 3: Normalized A
%
% See also amg, amgdoc
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Setup
if nargin == 0
    load lakemesh
    modification = 1; option = [];
end
if ~exist('modifcationn','var')
    modification = 1;
end
if ~exist('option','var')
    option = [];
end

%% Test
time = zeros(5,1); N = zeros(5,1); iter = zeros(5,1); err = zeros(5,1);
for k = 1:5
    [A,M] = assemblematrix(node,elem);
    switch modification
        case 0 % Neumann problem
            A = A(2:end,2:end);
        case 1 % Dirichlet problem
            [bdNode,bdEdge,isBdNode] = findboundary(elem); %#ok<*ASGLU>
            A = A(~isBdNode,~isBdNode);
        case 2 % Add mass matrix from time discretization
            h = 1/sqrt(size(elem,1)); %#ok<NASGU>
            dt = eval(varargin{1});
            A = A + M/dt;
        case 3             
            D = spdiags(1./sqrt(diag(A)),0,size(A,1),size(A,1));
            A = D*A*D;
    end
    N(k) = size(A,1);
    b = ones(N(k),1)/N(k);
    [x,info] = amg(A,b,option);
    time(k) = info.solverTime;
    iter(k) = info.itStep;
    err(k) = info.stopErr;
    if N(k) > 2e5  
        break
    end
    [node,elem] = uniformrefine(node,elem);
end

%% Output
time = time(1:k); N = N(1:k); iter = iter(1:k); err = err(1:k);
errHist = info.error;
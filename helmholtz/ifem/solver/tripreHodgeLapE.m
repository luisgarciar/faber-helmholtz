function [sigma,u,info] =  tripreHodgeLapE(M,G,C,f,g,node,elem,bdFlag,option)
%% TRIPREHODGELAPE: solve the Hodge Laplacian by a triangular preconditioner 
%
% [sigmau,relres,itStep] =  tripreHodgeLapE(M,B,C,f,g,node,elem,bdFlag,option)
% solves saddle point system discretized from mixed FEM
%
%      |M  G'| |sigma|  = |f|
%      |G -C | |u|      = |g|
%
% use MINRES method with triangular preconditioner where the Schur
% complement A = Gt*invM*G + C inverse is approximated by few V-cycles.
%
% Created by Jie Zhou on July 30, 2014.
%
% Rewritten by Long Chen on Feb 21, 2014.


t = cputime;
%% Parameters
if ~exist('option','var'), option = []; end
Ndof = length(f)+length(g); 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; tol = option.tol; maxIt = option.solvermaxit; 
printlevel = option.printlevel; 
d = size(node,2);

%% Set up auxiliary matrices
Nf = length(f); Ng = length(g);
if d == 2
    area = simplexvolume(node,elem);
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,...
                        [max(elem(:)),1]);
elseif d == 3
    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation presereved
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);    
end
DMinv = spdiags(1./Mvlump(option.isFreeNode),0,Nf,Nf);
% DMinv = spdiags(1./diag(M),0,N,N);
Gt = G';
Abar = G*DMinv*Gt + C;

%% Set up matrices for multigrid
setupOption.solver = 'NO';
setupOption.isFreeEdge = option.isFreeEdge;
[x,info,Ai,Bi,BBi,Res,Pro] = mgHodgeLapE(Abar,g,node,elem,bdFlag,setupOption); %#ok<*ASGLU>

%% Form a big matrix equation
bigA = [-M,Gt; G,C];
bigF = [f; g];

%% Preconditioned GMRES
% options for V-cycle of Schur complement
if strcmp(option.solver,'CG') % change default set up in mgoption 
   option.solver = 'Vcycle';
end
if isfield(option,'Vit')
   option.solvermaxIt = option.Vit;
else
   option.solvermaxIt = 1; 
end
option.setupflag = false;
option.x0 = zeros(Ng,1);
option.printlevel = 0;
% minres for the saddle point system
[x,flag,stopErr,itStep,err] = gmres(bigA,bigF,20,tol,maxIt,@tripre,[],x0);
itStep = (itStep(1)-1)*20 + itStep(2); % total iteration
sigma = x(1:Nf); 
u = x(Nf+1:end);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Triangular Preconditioned GMRES \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(bigA), option.solvermaxIt, itStep, stopErr, time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',stopErr);

%% Preconditioner 
    function e = tripre(r)
        r1 = r(1:Nf);
        r2 = r(Nf+1:end);
        e1 = -DMinv*r1;
        r2 = r2 - G*e1;
        e2 = mgHodgeLapE(Abar,r2,node,elem,bdFlag,option,Ai,Bi,BBi,Res,Pro);
        e  = [e1 + DMinv*(Gt*e2); e2];
    end
end
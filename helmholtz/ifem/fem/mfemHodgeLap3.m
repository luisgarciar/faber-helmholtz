function [err,time,solver,eqn] = mfemHodgeLap3(node,elem,pde,bdFlag,option,varargin)
%% MFEMHODGELAP3 solve Hodge Laplacian equation by various mixed finite element methods
%
%   MFEMHODGELAP3 computes approximations to the Hodge Laplacian equation
% 
% Example
%   HodgeLap3femrate
%
% See also mfemPoisson
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if ~exist('node','var') || ~exist('elem','var')
    [node,elem] = cubemesh([0,1,0,1,0,1],0.5);  % default mesh is a square
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = HodgeLaplacian3Edata1;                          % default data
end
if ~exist('bdFlag','var')
    bdFlag = setboundary3(node,elem,'Dirichlet'); 
end

%% Parameters
[elemType,maxIt,maxN,L0,refType,option] = mfemoption(option);

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
end

%% Initialize err
erruL2 = zeros(maxIt,1);  erruIuhH = zeros(maxIt,1); 
errsigmaL2 = zeros(maxIt,1); 
errsigmaH = zeros(maxIt,1); 
errsigmaIsigmah = zeros(maxIt,1); 
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1); h = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch elemType
        case 'ND0'  % lowest order Nedelec element
            [sigma,u,eqn,info] = HodgeLaplacian3E(node,elem,pde,bdFlag,option);
        case 'RT0'  % RT0 mixed FEM 
            [sigma,u,eqn,info] = HodgeLaplacian3F(node,elem,pde,bdFlag,option);
    end
    % compute error
    t = cputime;
    if isfield(pde,'sigma')
        errsigmaL2(k) = getL2error3(node,elem,pde.sigma,sigma);
        sigmaI = Lagrangeinterpolate(pde.sigma,node,elem);
        errSigma = sigma-sigmaI;
        errSigma = errSigma(eqn.freeNode);
        errsigmaIsigmah(k)=sqrt(errSigma'*eqn.Mv*errSigma);
    end
    if isfield(pde,'Dsigma')
        errsigmaH(k) = getH1error3(node,elem,pde.Dsigma,sigma);
    end
    if isfield(pde,'u')
        switch elemType
            case 'ND0'  % lowest order Nedelec element
                erruL2(k) = getL2error3NE(node,elem,pde.u,u);
                uI = edgeinterpolate(u,node,eqn.edge);
            case 'RT0'  % RT0 mixed FEM 
                erruL2(k) = getL2error3RT0(node,elem,pde.u,u);
                uI = faceinterpolate3(u,node,eqn.edge);
        end
        errU = u - uI;
        errU = errU(eqn.freeEdge);
        erruIuhH(k) = sqrt(errU'*eqn.C*errU);        
    end
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),erruL2(k), errsigmaL2(k));
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(u) + length(sigma);
    h(k) = 1./(sqrt(size(node,1))-1);
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
       figure(1);  
       showresult3(node,elem,sigma);    
    end
    if N(k) > maxN
        break;
    end
    % refine mesh
    t = cputime;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
    meshTime(k) = cputime - t;
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh2(h(1:k),erruIuhH(1:k),1,'-*','||D(u_I-u_h)||',...
               h(1:k),erruL2(1:k),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrateh3(h(1:k),errsigmaH(1:k),1,'-+','||D(\sigma - \sigma_h)||',... 
               h(1:k),errsigmaL2(1:k),1,'k-*','||\sigma - \sigma_h||',...
               h(1:k),errsigmaIsigmah(1:k),1,'m-+','||\sigma_I - \sigma_h||');
end

%% Output
err = struct('h',h(1:k),'N',N(1:k),'uL2',erruL2(1:k),'uIuhH',erruIuhH(1:k),...
             'sigmaL2',errsigmaL2(1:k),'sigmaH',errsigmaH(1:k),...
             'sigmaIsigmahL2',errsigmaIsigmah(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and CPU time
display('Table: Error')
colname = {'#Dof','h','||u-u_h||','||D(u_I-u_h)||','||sigma-sigma_h||','||D(sigma-sigma_h)||'};
displaytable(colname,err.N,[],err.h,'%0.2e',err.uL2,'%0.5e',err.uIuhH,'%0.5e',...
                     err.sigmaL2,'%0.5e',err.sigmaH,'%0.5e');
                 
display('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');   
              
display('Table: MG Iteration')
colname = {'#Dof','Iteration','Time'};
disptable(colname,solver.N,[],solver.itStep,'%2.0u',solver.time,'%4.2g');                 

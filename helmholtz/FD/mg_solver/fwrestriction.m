function R = fwrestriction(npf,dim,bc)
%% FWRESTRICTION Constructs the matrix corresponding to the full weight
%  restriction operator from a fine grid with npf interior points.
%  
%   Use:    R = fwrestriction(npf,dim,bc)  
%
%   Input:
%       npf:  number of interior points in 1-D fine grid (must be odd)
%       dim:  dimension (1 or 2)
%       bc:   'dir' for homogeneous Dirichlet boundary conditions
%             'som' for 1st order Sommerfeld boundary conditions 
%       
%   Output:
%       R:    restriction matrix, size npc x npf (1D Dirichlet)
%                                      npc^2 x npf^2 (2D Dirichlet)
%                                     (npc+2)^2 x (npf+2)^2 (2D Sommerfeld)
%   Author: Luis Garcia Ramos, 
%           Institut fur Mathematik, TU Berlin
%   Version 2.0, Sep 2016
%  Dirichlet and Sommerfeld boundary conditions
%  Works only on uniform grids: Modify this
%
%%
switch dim
    case 1
        npc = round(npf/2)-1; %number of interior points in coarse grid
        R = 2*lininterpol(npc,1,bc)';
    %   y = zeros(npf,1); y(1:3,1) = [1;2;1];
    %   R = gallery('circul',y');
    %   R = 0.25*sparse(R(1:2:(npf-2),:));
           
    case 2
        switch bc
            case 'dir'
                %npc = round(npf/2)-1;  
                y   = zeros(npf,1); y(1:3,1) = [1;2;1];
                R   = gallery('circul',y)';
                R   = 0.25*sparse(R(:,1:2:(npf-2))); %1D operator%
                R   = kron(R,R)';  %2D operator  
                
            case 'som'
                assert(mod(npf,2)==1,'number of interior points must be odd')
                npc = round((npf-1)/2);
                npff  = npf+2; % total number of points with endpoints
                npcc  = npc+2;
                R = 0.25*lininterpol(npc,dim,'som')';
                
                case 'som1'
                assert(mod(npf,2)==1,'number of interior points must be odd')
                npc = round((npf-1)/2);
                npff  = npf+2; % total number of points with endpoints
                npcc  = npc+2;
                R = 0.25*lininterpol(npc,dim,'som1')';
                 
        end      
      
end




         
%                 R     = sparse(npcc^2,npff^2);
% 
%                 %The restriction matrix is filled by rows 
%                 %(change this later!)
%                 %indc: row index (coarse grid)
%                 %indf: column index (fine grid)
% 
%                 %(0,0)- South-West Corner
%                 indc=1; indf=1; 
%                 R(indc,indf)=4; R(indc,indf+1)=4;
%                 R(indc,indf+npff)=4; R(indc,indf+npff+1)=4;
%                 
%                 %(1,0)-  South-East Corner
%                 indc=npcc; indf=npff;
%                 R(indc,indf)=4;  R(indc,indf-1)=4;
%                 R(indc,indf+npff)=4; R(indc,indf+npff-1)=4;
%                 
%                 %(0,1)- North-West Corner
%                 indc=npcc*(npcc-1)+1; indf=npff*(npff-1)+1;
%                 R(indc,indf)=4;   R(indc,indf-npf)=4;
%                 R(indc,indf+1)=4; R(indc,indf-npf+1)=4;
%                 
%                 %(1,1)- North East Corner
%                 indc=npcc^2; indf=npff^2;
%                 R(indc,indf)=4; R(indc,indf-1)=4;
%                 R(indc,indf-npff)=4; R(indc,indf-npff-1)=4;
%                              
%                 %South boundary y=0
%                 for indc=2:(npcc-1)
%                     indf=2*indc-1;
%                     R(indc,indf)=4; R(indc,indf+npff)=4;
%                     R(indc,indf+1)=2; R(indc,indf-1)=2;
%                     R(indc,indf+npff+1)=2; R(indc,indf+npff-1)=2; 
%                 end
%                 
%                 %East boundary x=1
%                 for i=2:(npcc-1)
%                     indc=i*npcc; indf=(2*i-1)*npff;
%                     R(indc,indf)=4; R(indc,indf-1)=4;
%                     R(indc,indf-1+npff)=2; R(indc,indf-1-npff)=2;
%                     R(indc,indf-npff)=2; R(indc,indf+npff)=2;
%                 end
%                 
%                 %North boundary y=1
%                 for i=2:(npcc-1)
%                    indc=(npcc-1)*npcc+i; indf=(npff-1)*npff+(2*i-1);
%                    R(indc,indf)=4; R(indc,indf-npf)=4; 
%                    R(indc,indf+1)=2; R(indc,indf-1)=2; 
%                    R(indc,indf-npff+1)=2; R(indc,indf-npff-1)=2;
%                 end
%                 
%                 %West boundary x=0
%                 for i=2:(npcc-1)
%                    indc=npcc*(i-1)+1; indf=npff*2*(i-1)+1;
%                    R(indc,indf)=4; R(indc,indf+1)=4;
%                    R(indc,indf+npff)=2;  R(indc,indf-npff)=2;
%                    R(indc,indf+1+npff)=2;  R(indc,indf+1-npff)=2;        
%                 end
%                              
%                 %Interior points
%                 for i=2:(npcc-1)
%                     for j=2:(npcc-1)
%                         indc=i+npcc*(j-1);
%                         ii=2*i-1; jj=2*j-1;
%                         indf=ii+npff*(jj-1);
%                         
%                         R(indc,indf)=4;
%                         R(indc,indf+npf)=2; R(indc,indf-npf)=2;
%                         R(indc,indf+1)=2; R(indc,indf-1)=2;
%                         R(indc,indf-1-npf)=1; R(indc,indf-1+npf)=1;
%                         R(indc,indf+1-npf)=1; R(indc,indf+1+npf)=1;
                     
end



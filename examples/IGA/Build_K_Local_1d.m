function [ Ke ] = Build_K_Local_1d( dR_dx,Jmod,E,A)
% [ Ke ] = Build_K_Local( dR_dx,Jmod,D,nen )
%-------------------------------------------------------------
% PURPOSE:
% Calculate the contribution from the current integration point 
% to the linear elastic stiffness matrix for a solid 
% NURBS element.
%
% INPUT: dR_dx = vector of basis function derivatives (nen x 3)
%
%        Jmod  = Jacobian determinant
%
%
%        nen   = number of local basis functions
%
% OUTPUT: Ke = element stiffness matrix contribution (nen*3 x nen*3)
%-------------------------------------------------------------

% Generate B matrix
B=dR_dx(:,1);

% Contrubution to element stiffness matrix
Ke = E*A*B*B'*Jmod;

end


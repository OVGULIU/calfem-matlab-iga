function [ Fe ] = Build_F_Local_1d( R,Jmod,force)
% [ Fe ] = Build_F_Local_1d( R,Jmod,force)
%-------------------------------------------------------------
% PURPOSE:
% Calculate the contribution from the current integration point 
% to the linear elastic residual vector for a solid 
% NURBS element.
%
% INPUT: R = vector of basis function  (nen x 1)
%
%        Jmod  = Jacobian determinant
%
%
%        nen   = number of local basis functions
%
% OUTPUT: Fe = element residual vector contribution (nen x 1)
%-------------------------------------------------------------

% Generate R matrix
R=R(:,1);

% Contrubution to element residual vector
Fe = force*R*Jmod;

end


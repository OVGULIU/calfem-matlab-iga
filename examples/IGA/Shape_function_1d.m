function [ R,dR_dx,J ] = Shape_function_1d( GP,e,deg,B,KV,INC,IEN)
% [ R,dR_dx,J ] = Shape_function( GP,e,deg,B,KV,INN,IEN)
%-------------------------------------------------------------
% PURPOSE:
% For a solid NURBS element, calculate the the vector of local 
% shape functions R and an array of their derivatives dR_dx at
% the current gauss point, and the Jacobian determinant J.
% 
% INPUT: The quadrature point GP, containing:
%        GP.Xi_tilde, GP.Eta_tilde, and GP.Zeta_tilde 
%
%        element number e, 
%
%        polynomial orders deg, containing:
%        deg.p, deg.q, and deg.r
%
%        control net B{n,m,l}(x,y,z,w), with a cell structure 
%        with n points in Xi direction, m points in Eta 
%        direction, and l points in Zeta direction. Weights 
%        are stored as the fourth component.
%
%        knot vectors KV, containing the clamped non-uniform 
%        knot vectors: 
%        KV.Xi, KV.Eta, and KV.Zeta
%
%        and connectivity arrays INC and IEN.
% 
% OUTPUT: R = vector of basis functions (nen x 1)
%         dR_dx = vector of basis function derivatives (nen x 3)
%         J = Jacaboian determinant (1)
%-------------------------------------------------------------

p = deg.p; 

% number of local basis functions:
nen = (p+1); 
                          
% NURBS coordinates; convention consistent with Algorithm 7
ni = INC(IEN(1,e),1);

% Calculate parametric coordinates from parent element coordinates
% Knot vectors KV.Xi, KV.Eta, and KV.Zeta and
% parent element coordinates xi_tilde, eta_tilde, zeta_tilde
% are given as input
xi = ((KV.Xi(ni+1)-KV.Xi(ni))*GP.xi_tilde ...
+ (KV.Xi(ni+1)+KV.Xi(ni))) / 2;


% Calculate univariate B-spline functions using (2.1) and (2.2)
% and their derivatives using (2.12)
N1 = Der1BasisFun(ni-1,xi,p,KV.Xi)'; % xi-dir.
N = N1(:,1);
dN_dxi = N1(:,2);

clear N1 

% Build numerators and denominators (in local numbering)
x = zeros(1,nen);
R = zeros(nen,1);   % Array of trivariate NURBS basis functions
dR_dxi = zeros(nen,3); % Trivariate NURBS function derivatives
                       % w.r.t. parametric coordinates
loc_num = 0; % Local basis function counter


        for i = 0 : p
            loc_num = loc_num + 1;
            
            R(loc_num) = N(p+1-i)* B{ni-i,1,1}(4); % Function numerator (N*M*L*w)
            
            % Get coordinates in local numbering
            x(loc_num) = B{ni-i,1,1}(1);
            %w(loc_num) = B{ni-i,nj-j,nk-k}(4);
            
            dR_dxi(loc_num,1) = dN_dxi(p+1-i)* B{ni-i,1,1}(4); % Derivative numerator (dN*M*L*w)
            

        end

        
W = sum(R); % Function denominator (Sum(N*M*L*w))
dW_dxi = sum(dR_dxi(:,1)); % Derivative denominator (Sum(dN*M*L*w))
            
% Divide by denominators to complete definitions of functions
% and derivatives w.r.t. parametric coordinates
dR_dxi = (dR_dxi*W - R*dW_dxi) / W^2;
R = R/W;

% Gradient of mapping from parameter space to physical space
dx_dxi = x* dR_dxi;

% Compute derivatives of basis functions
% with respect to physical coordinates
dR_dx = dR_dxi/dx_dxi; %dR_dxi * inv(dx_dxi)

% Gradient of mapping from parent element to parameter space
%dxi_dtildexi=zeros(3); % Derivative of parametric coordinates
                       % w.r.t. parent element coordinates
dxi_dtildexi = (KV.Xi(ni+1)-KV.Xi(ni))/2;

% Compute the Jacobian
J_mat = dx_dxi(1,1)*dxi_dtildexi;

% Compute Jacobian determinant
J = det(J_mat);

end


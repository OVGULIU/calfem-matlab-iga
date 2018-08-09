clear 
clc

%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Developed By:                                            %
%         %          Zhengkun Liu                          %
%         % email:   zhengkun.liu@ovgu.de                  %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set mesh [linear parameterization, continuity C^(p-k)]
k = 1; % default

% Polynomial degree and no. of elements
p = 6;
nel = 50; % summe of nodes is p+nel

% Length of the bar
lang = 10;

% p-h-k refinement using NUBRS Toolbox
knot      = [zeros(1,p+1),ones(1,p+1)];
h         = 1/nel;
ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
[cp,knot] = bspkntins(p,0:1/p:1,knot,ins);
cp = cp*lang;
ncp = length(cp);

% Change of the Format for CALFEM-IgA
CP = zeros(ncp,4);

for i=1:ncp
   CP(i,1)=cp(i);
   CP(i,4)=1;
end



B = cell(ncp,1,1);

for i=1:ncp
    B{i,1,1}=CP(i,:);
end

% No. of Gauss Points
ngauss=p+1;

Xi = knot';

k_ = length(Xi);

%%%%%%%%

% Number of weights
n = size(B,1);
 
% Order of basis
p = k_-n-1;
 
deg.p = p;  clear p 

% 
 
%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; clear Xi 
 
% Build connectivity arrays
nel = n-deg.p; % number of elements
nnp = n; % number of global basis functions
nen = deg.p+1; % number of local basis functions
ndof = nnp; % number of global degrees of freedom
ldof = nen; % number of local degrees of freedom
% Build connectivity arrays
[INN,IEN] = BldINCIEN( deg.p,n ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
ID = reshape(1:max(max(IEN)),1,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed)
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),nen,1);
end

%% Material parameters:
E_Y = 1e5;
A=1; 
% external loading
Force = 500;
%% Gauss-Legendre quadrature points:
[ gp_x,w_x ] = getGP( deg.p );
NQUADx = size(gp_x,2);
 
%% Stiffness matrix and load vector computation
 
% Element loop
K = zeros(ndof); % Needs to be changed to sparse for large problems!!

F = zeros(ndof,1);
for e = 1 : nel
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    
    % Check if element has zero measure
    if KV.Xi(ni+1) == KV.Xi(ni)
        continue
    end
    
    Ke = zeros(nen);
    Fe = zeros(nen,1);
    for i = 1 : NQUADx % Loop trough Gauss points
                % Gauss point
                GP.xi_tilde = gp_x(i);
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function_1d( GP,e,deg,B,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i);
                
                % Build Ke
                [ Ke_ ] = Build_K_Local_1d( dR_dx,Jmod,E_Y,A);
                [ Fe_ ] = Build_F_Local_1d( R,Jmod,Force);
                Ke = Ke + Ke_;
                Fe = Fe + Fe_;
    end
 
    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
    F(idx) = F(idx)+Fe;
end
 

%   Calfem FE-solver 
      
     [nd,nd]=size(K);
      fdof=[1:nd]';
%
     d=zeros(size(fdof));
     Q=zeros(size(fdof));
%
     bc = [1 0; nd 0]; 
     pdof=bc(:,1);
     dp=bc(:,2);
     fdof(pdof)=[];
%
     s=K(fdof,fdof)\(F(fdof)-K(fdof,pdof)*dp);
%
     d(pdof)=dp;
     d(fdof)=s;




% post-process 
plot(CP(:,1),d,'--ok')




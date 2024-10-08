function [D] = GetStressStag(D,N)
% ======================================================================= %
% Calculate the strain rate tensor for a staggered grid velocity grid     %
% ----------------------------------------------------------------------- %
%                                                                         %
%       .----o----.----o----.                                             %
%       |         |         |       . - eta_s, e_{xz,zx}, t_{xz,zx}       %
%       +    *    +    *    +       * - eta_n, e_{xx,zz}, t_{xx,zz}, P    %
%       |         |         |       o - v_z                               %
%       .----o----.----o----.       + - v_x                               %
%       |         |         |                                             %
%       +    *    +    *    +       e_{i,j} = 0.5*(dv_i/dx_j + dv_j/dx_i) %
%       |         |         |       t_{i,j = 2 * eta * e_{i,j}            %
%       .----o----.----o----.                                             %
%                                                                         %
% ----------------------------------------------------------------------- %
% LF - 22.01.24 -                                                         %
% ======================================================================= %

% Define Field Size ===================================================== %
D.etac      =   zeros(N.nz-1,N.nx-1);
D.tauxx     =   zeros(N.nz-1,N.nx-1);
D.tauzz     =   zeros(N.nz-1,N.nx-1);
D.tauxz     =   zeros(N.nz,N.nx);
D.tauII     =   zeros(N.nz,N.nx);
% ======================================================================= %

% Calculate Viscosity on Centroids ====================================== %
A1          =   zeros(N.nz,N.nx);

A1(1:N.nz-1,:)    =   ...
    (D.eta(2:N.nz,:) + D.eta(1:N.nz-1,:))/2;

D.etac      =   (A1(1:N.nz-1,2:N.nx) + A1(1:N.nz-1,1:N.nx-1)) /2;
% ======================================================================= %

% Stress Tensor ========================================================= %
% Nodes -----------------------------------------
D.tauxz         =   2 .* D.eta .* D.exz;
% Centroids -------------------------------------
D.tauxx         =   2 .* D.etac .* D.dvxdx;
D.tauzz         =   2 .* D.etac .* D.dvzdz;
% ======================================================================= %

% 2nd Invariant on Nodes ================================================ %
D.tauII     =   sqrt(0.5.*(c2n(D.tauxx).^2 + c2n(D.tauzz).^2 + ...
    2.*D.tauxz.^2));
% ======================================================================= %
end

function A2 = c2n(A0)

A1      = zeros(size(A0,1)+1,size(A0,2));

A1(:,:) = [1.5*A0(1,:)-0.5*A0(2,:); 
            (A0(2:end,:)+A0(1:end-1,:))/2; 
            1.5*A0(end,:)-0.5*A0(end-1,:)];
        
A2      = zeros(size(A1,1),size(A1,2)+1);
A2(:,:) = [1.5*A1(:,1)-0.5*A1(:,2), ...
    (A1(:,2:end)+A1(:,1:end-1))/2, ...
    1.5*A1(:,end)-0.5*A1(:,end-1)];
end






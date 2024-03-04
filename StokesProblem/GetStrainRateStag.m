function [D] = GetStrainRateStag(D,N)
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
D.dvxdx     =   zeros(N.nz-1,N.nx-1);
D.dvzdz     =   zeros(N.nz-1,N.nx-1);
D.dvxdz     =   zeros(N.nz,N.nx);
D.dvzdx     =   zeros(N.nz,N.nx);
% ======================================================================= %

% Calculate Velocity and Strain Rate Tensor ============================= %
% On the corresponding positions as shown above.                          %
% Velocity -------------------------------------------------------------- %
% Centroids -----------------
D.dvxdx  =   (D.vx(1:N.nz-1,2:N.nx) - D.vx(1:N.nz-1,1:N.nx-1))/N.dx;
D.dvzdz  =   (D.vz(2:N.nz,1:N.nx-1) - D.vz(1:N.nz-1,1:N.nx-1))/N.dz;
% Nodes ---------------------
D.dvxdz(2:N.nz-1,1:N.nx)    =   ...
    (D.vx(2:N.nz-1,1:N.nx) - D.vx(1:N.nz-2,1:N.nx))/N.dz;
% Linear extrapolation ------
D.dvxdz(1,:)            =   D.dvxdz(2,:) + (D.dvxdz(3,:) - D.dvxdz(2,:));
D.dvxdz(N.nz,:)         =   ...
    D.dvxdz(N.nz-1,:) + (D.dvxdz(N.nz-2,:) - D.dvxdz(N.nz-1,:));

D.dvzdx(1:N.nz,2:N.nx-1)    =   ...
    (D.vz(1:N.nz,2:N.nx-1) - D.vz(1:N.nz,1:N.nx-2))/N.dx;
% Linear extrapolation ------
D.dvzdx(:,1)            =   D.dvzdx(:,2) + (D.dvzdx(:,3) - D.dvzdx(:,2));
D.dvzdx(:,N.nx)         =   ...
    D.dvzdx(:,N.nx-1) + (D.dvzdx(:,N.nx-2) - D.dvzdx(:,N.nx-1));

% Strain Rate ----------------------------------------------------------- %
D.exx       =   D.dvxdx;
D.ezz       =   D.dvzdz;
D.exz       =   0.5*(D.dvxdz + D.dvzdx);
% ======================================================================= %

% Invariant on the nodes ================================================ %
D.eII   =   sqrt( 0.5 .* (c2n(D.exx).^2 + c2n(D.ezz).^2 + 2.*D.exz.^2));
% ======================================================================= %

end

function A2 = c2n(A0)

A1      = zeros(size(A0,1)+1,size(A0,2));

A1(:,:) = [1.5*A0(1,:)-0.5*A0(2,:); 
            (A0(2:end,:)+A0(1:end-1,:))/2; 
            1.5*A0(end,:)-0.5*A0(end-1,:)];
        
A2      = zeros(size(A1,1),size(A1,2)+1);
A2(:,:) = [1.5*A1(:,1)-0.5*A1(:,2), (A1(:,2:end)+A1(:,1:end-1))/2, 1.5*A1(:,end)-0.5*A1(:,end-1)];
end





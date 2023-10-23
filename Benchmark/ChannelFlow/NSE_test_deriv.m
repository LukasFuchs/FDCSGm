%% ====================================================================== %
% Program to calcutlate the first derivative of the function              %
% df/dx = d/dx eta * dv/dx, i.e., the 1-D stokes equation for variable    % 
% viscosity and a horizontal pressure gradient using different finite     % 
% difference schemes. The schemes used hear are:                          %
%   -  
% 
% ======================================================================= %
% LF - 10/19/23 -                                                         %
% ======================================================================= %
clear
clc
clf
%% Scaling constants ==================================================== %
Hsc         =   400e3;              % Reference height [ m ]
vsc         =   1.58e-9;            % Reference velocity [ m/s ]
etasc       =   1e20;               % Reference viscosity [ Pa s ]
% ======================================================================= %
%% Model constants ====================================================== %
nz          =   51;                 % Number of grid nodes
H           =   400e3;              % Model height [ m ]
z0          =   0;                  % Surface coordinate [ m ]
v0          =   1.58e-9;            % Velocity at the top [ m/s ]
v1          =   0;                  % Velocity at the bottom [ m/s ]
eta0        =   20;                 % log10(eta) at the bottom
eta1        =   25;                 % log10(eta) at the top
m           =   10^eta1/10^eta0;    % Viscosity ration top/bottom
k           =   logspace(eta0,eta1,nz)';
dPdx        =   -0.0;               % non dimensional pressure gradient
% ======================================================================= %
%% Scaling parameters =================================================== %
H           =   H/Hsc; 
v0          =   v0/vsc;
dz          =   H/(nz-1);           % Grid spacing
z           =   (z0:dz:H)';         % Depth 
k           =   k./etasc; 
eta         =   10^eta0.*...        % Logarithmic Viscosity profile 
                    exp(log(m).*z./H);
eta         =   eta./etasc;
dfdx        =   ones(nz,1).*dPdx; 
dfdx(1)     =   v1; 
dfdx(nz)    =   v0; 
% ======================================================================= %
%% Analytical shear velocity ============================================ %
if m == 1
    v_ana   =   -dPdx/2/eta(1)*(H.*z - z.^2) +...
                    v0.*z./H;
else
    v_ana   =   -dPdx*H/(eta(1))/log(m).*...
                    (m.^(-z./H)/(m-1).*(z.*(m-1) + H) - H/(m-1)) - ...
                    m.^(-z./H).*m.*v0./(m-1) + ...
                    v0*m/(m-1);
end
% ======================================================================= %
%% staggered grid ======================================================= %
a1      =   zeros(nz,1);
b1      =   zeros(nz,1);
c1      =   zeros(nz,1); 
for i = 2:(nz-1)
    a1(i-1) = (eta(i-1) + eta(i))/2/dz^2;
    b1(i)   = -(eta(i-1) + 2*eta(i) + eta(i+1))/2/dz^2;
    c1(i+1) = (eta(i) + eta(i+1))/2/dz^2;
end
b1(1)   =   1; 
b1(nz)  =   1; 
abc1    =   [a1,b1,c1];
A1      =   spdiags(abc1,[-1 0 1],nz,nz);
v_num   =   A1\dfdx;
% STRESS CALCULATION USING MATRIX NOTATION ------------------------------ %
a12     =   zeros(nz,1); 
b12     =   zeros(nz,1); 
c12     =   zeros(nz,1); 
for i = 2:(nz-1)
    a12(i-1)    =   -(eta(i-1) + eta(i))/4/dz;
    b12(i)      =   (eta(i-1) - eta(i+1))/4/dz;
    c12(i+1)    =   (eta(i) + eta(i+1))/4/dz;
end
b12(1)      =   -(eta(1) + eta(2))/2/dz;
b12(end)    =   (eta(end) + eta(end-1))/2/dz;
c12(2)      =   (eta(1) + eta(2))/2/dz;
a12(end-1)  =   -(eta(end) + eta(end-1))/2/dz;
abc12       =   [a12,b12,c12];
A12         =   spdiags(abc12,[-1 0 1],nz,nz);
tauxz1      =   A12*v_num;
% ======================================================================= %
%% Chain rule, forward differences  ====================================== %
a2      =   zeros(nz,1);
b2      =   zeros(nz,1);
c2      =   zeros(nz,1);
for i = 2:(nz-1)
    a2(i-1) =   (eta(i-1) + 4*eta(i) - eta(i+1))/4/dz^2;
    b2(i)   =   -2*eta(i)/dz^2;
    c2(i+1) =   -((eta(i-1) - 4*eta(i) - eta(i+1))/4/dz^2);
end
b2(1)   =   1;
b2(nz)  =   1;
abc2    =   [a2,b2,c2];
A2      =   spdiags(abc2,[-1 0 1],nz,nz);
v_num2  =   A2\dfdx;
% % STRESS CALCLUATIONS USING MATRIX NOTATION ============================= %
% a22     =   zeros(nz,1);
% b22     =   zeros(nz,1);
% for i=2:(nz-1)
%     a22(i-1)    =   eta(i)/2/dz;
%     b22(i+1)    =   eta(i)/2/dz;
% end
% abc22   =   [a22,b22];
% A22     =   spdiags(abc22,[-1 1],nz,nz);
% tauxz2  =   A22*v_num2;
% ======================================================================= %
%% Centered without staggered grid points =============================== %
a3      =   zeros(nz,1);    b3  =   zeros(nz,1);    c3  =   zeros(nz,1);
as      =   zeros(nz,1);    bs  =   zeros(nz,1); 
a3(nz-3)=   eta(nz-2)/4/dz^2;
as(1)   =   eta(1)/2/dz^2;
b3(1)   =   1; 
b3(2)   =   -(eta(3)+2*eta(2))/4/dz^2;
b3(nz-1)=   -(eta(nz-1)+2*eta(nz))/4/dz^2;
b3(nz)  =   1; 
bs(nz)  =   eta(nz)/2/dz^2;
c3(4)   =   eta(3)/4/dz^2; 
for i = 3:(nz-2)
    a3(i-2) =   eta(i-1)/4/dz^2;
    b3(i)   =   -(eta(i-1)+eta(i+1))/4/dz^2;
    c3(i+2) =   eta(i+1)/4/dz^2;
end
abc3        =   [a3,as,b3,bs,c3];
A3          =   spdiags(abc3,[-2 -1 0 1 2],nz,nz);
v_num3      =   A3\dfdx;
% STRESS CALCULATION USING MATRIX NOTATION ------------------------------ %
a32         =   zeros(nz,1);    b32         =   zeros(nz,1); 
c32         =   zeros(nz,1);    d32         =   zeros(nz,1); 
e32         =   zeros(nz,1); 
for i = 2:(nz-1)
    a32(i-1) =  -eta(i)/2/dz;
    c32(i+1) =  eta(i)/2/dz;
end
b32(1)      =   -3*eta(1)/2/dz;
c32(2)      =   2*eta(1)/dz;
d32(3)      =   -eta(1)/2/dz;
b32(end)    =   3*eta(end)/2/dz;
a32(end-1)  =   -2*eta(end)/dz;
e32(end-2)  =   eta(end)/2/dz;
abc32       =   [e32,a32,b32,c32,d32];
A32         =   spdiags(abc32,[-2 -1 0 1 2],nz,nz);
tauxz3      =   A32*v_num3;
% ======================================================================= %
%% Fourth order, centered without staggered grid points ================= %
a4      =   zeros(nz,1);    b4  =   zeros(nz,1);    c4  =   zeros(nz,1);
d4      =   zeros(nz,1);    e4  =   zeros(nz,1);    f4  =   zeros(nz,1);
g4      =   zeros(nz,1);    h4  =   zeros(nz,1);    j4  =   zeros(nz,1);
m4      =   zeros(nz,1);    z4  =   zeros(nz,1);
e4(1)   =   1;
e4(nz)  =   1;
a       =   1/12;
b       =   -2/3;
ap      =   -25/12;
bp      =   4;
cp      =   -3;
dp      =   4/3;
ep      =   -1/4;
a2p     =   -1/4;
b2p     =   -5/6;
c2p     =   3/2;
d2p     =   -1/2;
e2p     =   1/12;
for i = 5:(nz-4)
    a4(i-4) =   a^2*eta(i-2);
    b4(i-3) =   a*b*(eta(i-2) + eta(i-1));
    c4(i-2) =   b^2*eta(i-1);
    d4(i-1) =   -a*b*(eta(i-2) + eta(i+1));
    e4(i)   =   -(a^2*(eta(i-2) + eta(i+2)) + b^2*(eta(i-1) + eta(i+1)));
    f4(i+1) =   -a*b*(eta(i-1) + eta(i+2));
    g4(i+2) =   b^2*eta(i+1);
    h4(i+3) =   a*b*(eta(i+1) + eta(i+2));
    j4(i+4) =   a^2*eta(i+2);
end
b4(1)       =   a*(a2p*eta(2) + b*eta(3));
c4(2)       =   a*b2p*eta(2) + b^2*eta(3);
d4(3)       =   a*(eta(2)*c2p - b*eta(5));
e4(4)       =   a*(eta(2)*d2p - a*eta(6)) - b^2*(eta(3) + eta(5));
f4(5)       =   a*(eta(2)*e2p - b*(eta(3) + eta(6)));
g4(6)       =   b^2*eta(5);
h4(7)       =   a*b*(k(5) + eta(6));
j4(8)       =   a^2*eta(6);
a4(nz-7)    =   a^2*eta(nz-5);
b4(nz-6)    =   a*b*(eta(nz-5) + eta(nz-4));
c4(nz-5)    =   b^2*eta(nz-4);
d4(nz-4)    =   -a*(b*(eta(nz-5) + eta(nz-2)) - e2p*eta(nz-1));
e4(nz-3)    =   -(a*(a*eta(nz-5) - d2p*eta(nz-1)) + b^2*(eta(nz-4) + eta(nz-2)));
f4(nz-2)    =   a*(c2p*eta(nz-1) - b*eta(nz-4));
g4(nz-1)    =   b^2*eta(nz-2) + a*b2p*eta(nz-1);
h4(nz)      =   a*(b*eta(nz-2) + a2p*eta(nz-1));
c4(1)       =   a*ap*eta(1) + b*a2p*eta(2);
d4(2)       =   a*(bp*eta(1) - b*eta(4)) + b*b2p*eta(2);
e4(3)       =   a*(cp*eta(1) - a*eta(5)) + b*(c2p*eta(2) - b*eta(4));
f4(4)       =   a*(dp*eta(1) - b*eta(5)) + b*d2p*eta(2);
g4(5)       =   a*ep*eta(1) + b*(e2p*eta(2) + b*eta(4));
h4(6)       =   a*b*(eta(4) + eta(5));
j4(7)       =   a^2*eta(5);
a4(nz-6)    =  a^2*eta(nz-4);
b4(nz-5)    =  a*b*(eta(nz-4) + eta(nz-3));
c4(nz-4)    =  b*(b*eta(nz-3) + e2p*eta(nz-1)) + a*ep*eta(nz);
d4(nz-3)    =  -(b*(a*eta(nz-4) - d2p*eta(nz-1)) - a*dp*eta(nz));
e4(nz-2)    =  -(a*(a*eta(nz-4) - cp*eta(nz)) + b*(b*eta(nz-3) - c2p*eta(nz-1)));
f4(nz-1)    =  -(a*(b*eta(nz-3) - bp*eta(nz)) - b*b2p*eta(nz-1));
g4(nz)      =  b*a2p*eta(nz-1) + a*ap*eta(nz);
d4(1)       =   a2p*(ap*eta(1) + b2p*eta(2)) + a*c2p*eta(3);
e4(2)       =   bp*a2p*eta(1) + b2p^2*eta(2) + b*c2p*eta(3) + a*d2p*eta(4) ;
f4(3)       =   cp*a2p*eta(1) + b2p*c2p*eta(2) + b*d2p*eta(4) + a*e2p*eta(5);
g4(4)       =   dp*a2p*eta(1) + b2p*d2p*eta(2) - b*(c2p*eta(3) - e2p*eta(5));
h4(5)       =   ep*a2p*eta(1) + b2p*e2p*eta(2) - a*c2p*eta(3) - b*d2p*eta(4);
j4(6)       =   -(a*d2p*eta(4) + b*e2p*eta(5));
m4(7)       =   -a*e2p*eta(5);
z4(nz-6)    =   -(a*e2p*eta(nz-4));
a4(nz-5)    =   -(b*e2p*eta(nz-4) + a*d2p*eta(nz-3));
b4(nz-4)    =   -(b*d2p*eta(nz-3) + a*c2p*eta(nz-2) - b2p*e2p*eta(nz-1) -...
                    ep*a2p*eta(nz));
c4(nz-3)    =   b*(e2p*eta(nz-4) - c2p*eta(nz-2)) + b2p*d2p*eta(nz-1) + ...
                    dp*a2p*eta(nz);
d4(nz-2)    =   a*e2p*eta(nz-4) + b*d2p*eta(nz-3) + c2p*b2p*eta(nz-1) + ...
                    cp*a2p*eta(nz);
e4(nz-1)    =   a*d2p*eta(nz-3) + b*c2p*eta(nz-2) + b2p^2*eta(nz-1) + ...
                    bp*a2p*eta(nz);
f4(nz)      =   a*c2p*eta(nz-2) + a2p*b2p*eta(nz-1) + ap*a2p*eta(nz);
z4          =   z4./dz^2;
a4          =   a4./dz^2;
b4          =   b4./dz^2;
c4          =   c4./dz^2;
d4          =   d4./dz^2;
e4          =   e4./dz^2;
f4          =   f4./dz^2;
g4          =   g4./dz^2;
h4          =   h4./dz^2;
j4          =   j4./dz^2;
m4          =   m4./dz^2;
e4(1)       =   1; 
e4(nz)      =   1;
abc4        =   [z4,a4,b4,c4,d4,e4,f4,g4,h4,j4,m4];
A4          =   spdiags(abc4,[-5 -4 -3 -2 -1 0 1 2 3 4 5],nz,nz);
v_num4      =   A4\dfdx;
% STRESS CALCULATION USING MATRIX NOTATION ------------------------------ %
x42     =   zeros(nz,1);    z42     =   zeros(nz,1);    a42     =   zeros(nz,1);
b42     =   zeros(nz,1);    c42     =   zeros(nz,1);    d42     =   zeros(nz,1);
e42     =   zeros(nz,1);    f42     =   zeros(nz,1);    g42     =   zeros(nz,1);
for i = 3:(nz-2)
    a42(i-2)    =   eta(i)/12/dz;
    b42(i-1)    =   -2*eta(i)/3/dz;
    d42(i+1)    =   -b42(i-1);
    e42(i+2)    =   -a42(i-2);
end
for i = 2
    b42(i-1)    =   -eta(i)/4/dz;
    c42(i)      =   -5*eta(i)/6/dz;
    d42(i+1)    =   3*eta(i)/2/dz;
    e42(i+2)    =   -eta(i)/2/dz;
    f42(i+3)    =   eta(i)/12/dz;
end
for i = (nz-1)
   d42(i+1)     =   eta(i)/4/dz;
   c42(i)       =   5*eta(i)/6/dz;
   b42(i-1)     =   -3*eta(i)/2/dz;
   a42(i-2)     =   eta(i)/2/dz;
   z42(i-3)     =   -eta(i)/12/dz; 
end
for i = 1
    c42(i)      =   -25*eta(i)/12/dz;
    d42(i+1)    =   4*eta(i)/dz;
    e42(i+2)    =   -3*eta(i)/dz;
    f42(i+3)    =   4*eta(i)/3/dz;
    g42(i+4)    =   -eta(i)/4/dz;
end
for i = nz
    c42(i)      =   25*eta(i)/12/dz;
    b42(i-1)    =   -4*eta(i)/dz;
    a42(i-2)    =   3*eta(i)/dz;
    z42(i-3)    =   -4*eta(i)/3/dz;
    x42(i-4)    =   eta(i)/4/dz;
end
abc42   =   [x42,z42,a42,b42,c42,d42,e42,f42,g42];
A42     =   spdiags(abc42,[-4 -3 -2 -1 0 1 2 3 4],nz,nz); 
tauxz4  =   A42*v_num4;
% ======================================================================= %
%% Deviation from the analytical solution =============================== %
d_num1  =   sqrt((v_num-v_ana).^2./max(v_ana).^2);
d_num2  =   sqrt((v_num2-v_ana).^2./max(v_ana).^2);
d_num3  =   sqrt((v_num3-v_ana).^2./max(v_ana).^2);
d_num4  =   sqrt((v_num4-v_ana).^2./max(v_ana).^2);
% ======================================================================= %
%% Plot data ============================================================ %
nmsize  =   5; 
figure(1)
subplot(3,2,1)
plot(v_ana,z,'k-',...
    v_num,z,'k*',...
    v_num2,z,'rd',...
    v_num3,z,'gs',...
    v_num4,z,'mo','MarkerSize',nmsize)
legend('v_{ana}','v_{steg}','v_{Chain}','v_{centr2nd}','v_{centr4th}',...
    'Location','NorthWest')
ylabel('z','Interpreter','latex'); xlabel('v','Interpreter','latex')
title('Velocity','Interpreter','latex')
set(gca,'FontWeight','Bold',...
    'LineWidth',2,'FontSize',12,'TickLabelInterpreter','latex')
subplot(3,2,3) 
semilogx(eta,z,'LineWidth',2)
ylabel('z','Interpreter','latex')
xlabel('$$\eta$$','Interpreter','latex')
title('Viscosity','Interpreter','latex')
set(gca,'FontWeight','Bold',...
    'LineWidth',2,'FontSize',12,'TickLabelInterpreter','latex')
subplot(1,2,2)
semilogy(z,d_num1,'k*',...
    z,d_num2,'rd',...
    z,d_num3,'gs',...
    z,d_num4,'mo','MarkerSize',nmsize)
legend('$$\Delta_{stag}$$',...
        '$$\Delta_{Chain}$$',...
        '$$\Delta_{centr2nd}$$',...
        '$$\Delta_{centr4th}$$',...
        'Interpreter','latex','Location','Best')
xlabel('z','Interpreter','latex')
ylabel('$$\Delta$$','Interpreter','latex')
title('Square root error','Interpreter','latex')
set(gca,'FontWeight','Bold',...
    'LineWidth',2,'FontSize',12,'TickLabelInterpreter','latex')
subplot(3,2,5)
plot(tauxz1,z,'k',...           %tauxz2,z,'y--',...        %tauxz3,z,'r--',...
        tauxz4,z,'m--','LineWidth',2)
xlabel('$$\tau_{xz}$$','Interpreter','latex')
ylabel('z','Interpreter','latex')
title('Shear stress','Interpreter','latex')
% legend('stag','centr4th','Location','Best')
set(gca,'FontWeight','Bold',...
    'LineWidth',2,'FontSize',12,'TickLabelInterpreter','latex')

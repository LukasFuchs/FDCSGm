function [Py,D,ID,M,N,A] = SetUpFields(Py,N,M,T)
%% ------------------- Berechnung Numerischer Parameter ----------------- %
M.H         =   M.H*1000;       %   jetzt in [m]
M.L         =   -M.xmax*M.H;    %   Model Laenge  - [m]

N.nz1       =   (N.nz-1);
N.nx1       =   N.nx-1;
N.dz        =   M.H/(N.nz1);
N.dx        =   M.L/(N.nx1);

% Koordinaten des regulaeren Gitters
M.x         =   (0:N.dx:M.L);
M.z         =   (0:N.dz:M.H);
[M.X,M.Z]   =   meshgrid(M.x,M.z);

% Koordinaten des Zentrums eines FD-Elements
M.x1        =   (N.dx/2:N.dx:(M.L-N.dx/2));
M.z1        =   (N.dz/2:N.dz:(M.H-N.dz/2));
[M.X1,M.Z1] =   meshgrid(M.x1,M.z1);
A           =   [];
clear M.x M.z M.x1 M.z1
% ----------------------------------------------------------------------- %

% Physikalische Parameter auf dem regularen Gitter
D.rho      =   zeros(N.nz,N.nx);
D.C        =   zeros(N.nz,N.nx);
D.T        =   zeros(N.nz,N.nx);
D.Q        =   ones(N.nz,N.nx).*Py.Q0;

% "Staggered" Gitterpunkte
D.P         =   zeros(N.nz,N.nx);
D.vx        =   zeros(N.nz,N.nx);
D.vz        =   zeros(N.nz,N.nx);

D.Pi        =   zeros(N.nz,N.nx);
D.vxi       =   zeros(N.nz,N.nx);
D.vzi       =   zeros(N.nz,N.nx);

D.meanT     =   zeros(N.nz,T.itmax);
D.meanV     =   zeros(T.itmax,1);
D.Nus       =   zeros(T.itmax,1); 

% Interpolierte Parameter auf dem regulaeren Gitter
ID.vx       =   zeros(N.nz,N.nx);
ID.vz       =   zeros(N.nz,N.nx);
ID.v        =   zeros(N.nz,N.nx);
ID.vxo      =   [];
ID.vzo      =   [];
ID.indi     =   2:N.nz1;
ID.indj     =   2:N.nx1;

end
clear
clc

%% Model Geometry ======================================================= %
H       =   -400e3;             % [ m ]
vx0     =   5;                  % [ cm/a ]
vx0     =   vx0/100/31536000;   % [ m/s ]
eta     =   1e21;               % [ Pa s ]
dPdx    =   -1e2;               % [ Pa/m ]; z.B. -1e2 (Warum negativ?)

%% Numerische Parameter ================================================= %
nz      =   101; 
dz      =   H/(nz-1); 

% Koordinaten
z       =   H:-dz:0;
% Index der internen Gitterpunkte
indi    =   2:nz-1; 

%% Erstellung der Koeffizientenmatrix A 
diag            =   zeros(nz,3); 

% Bestimmung der Diagonalen 
diag(indi-1,1)  =   1/dz^2; 
diag(indi,2)    =   -2/dz^2; 
diag(indi+1,3)  =   1/dz^2; 

% Randbedingungen - no slip, d.h. konstante Geschwindigkeit
diag(1,2)       =   1; 
diag(nz,2)      =   1; 

A               =   spdiags(diag,[-1 0 1],nz,nz);

%% Definition der rechten Seite des Gleichungssystems
rhs             =   zeros(nz,1); 

rhs(indi)       =   dPdx/eta;
rhs(1)          =   0; 
rhs(nz)         =   vx0; 

vx              =   A\rhs; 

%% Analytische Loesung ================================================== %
vx_ana          =   1/2/eta*dPdx.*(z.^2 - H.*z) - vx0.*z./H + vx0;

%% Darstellung der Daten
figure(1)
clf
plot(vx,z./1e3,'LineWidth',2)
hold on
plot(vx_ana,z./1e3,'r--','LineWidth',2)
xlabel('v_x [ m/s ]'); ylabel('z [km]'); title('Geschwindigkeitsprofil')
set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2);

%% ===================================================================== %%
% ============================== ENDE =================================== %
% ======================================================================= %













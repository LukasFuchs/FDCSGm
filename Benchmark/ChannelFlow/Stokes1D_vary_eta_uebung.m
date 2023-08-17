clear
clc

%% Model Geometry ======================================================= %
H       =   -400e3;             % Maechtigkeit der Schicht [ m ]
vx0     =   2;                  % Geschwindigkeit oben [ cm/a ]
vx0     =   vx0/100/31536000;   % [ m/s ]
dPdx    =   -1e2;               % Druckgradient [ Pa/m ]; z.B. -7e1 (Warum negativ?)

%% Numerische Parameter ================================================= %
nz      =   ?;
dz      =   ?;

% Koordinaten
z       =   H:-dz:0;
% Index der internen Gitterpunkte
indi    =   2:nz-1;

%% Definition der Vikosität ============================================= %
eta0        =   22;                 % Power der Viskosität oben, e.g. 22
eta1        =   21;                 % Power der Viskosität unten, e.g. 21
m           =   10^eta1/10^eta0;    % Viskositätsverhältnis
eta         =   10.^eta0*exp(log(m).*z'./H);

%% Analytische Loesung ================================================== %
vx_ana      =   -dPdx*H/10^eta0/log(m)/(m-1).*...
    (z.*(m.^((-z+H)./H) - m.^(-z./H)) + H.*(m.^(-z./H) - 1)) + ...
    vx0/(m-1).*(m.^((-z+H)./H) - 1);

%% Erstellung der Koeffizientenmatrix A
diag            =   zeros(nz,3);

% Bestimmung der Diagonalen
diag(indi-1,1)  =   ?;
diag(indi,2)    =   ?;
diag(indi+1,3)  =   ?;

% Randbedingungen - no slip, d.h. konstante Geschwindigkeit
diag(1,2)       =   ?;
diag(nz,2)      =   ?;

A               =   spdiags(diag,[-1 0 1],nz,nz);

%% Definition der rechten Seite des Gleichungssystems
rhs             =   zeros(nz,1);

rhs(indi)       =   ?;
rhs(1)          =   ?;          % unten
rhs(nz)         =   ?;          % oben

vx              =   A\rhs;

%% Darstellung der Daten
figure(1)
clf
subplot(1,2,1)
plot(vx,z./1e3,'LineWidth',2)
hold on
plot(vx_ana,z./1e3,'r--','LineWidth',2)
xlabel('v_x [ m/s ]'); ylabel('z [km]'); title('Geschwindigkeitsprofil')
set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2);
subplot(1,2,2)
plot(eta,z./1e3,'LineWidth',2)
xlabel('\eta[ Pa s ]'); ylabel('z [km]'); title('Viskosity')
set(gca,'xscale','log','FontWeight','Bold','FontSize',15,'LineWidth',2);

%% ===================================================================== %%
% ============================== ENDE =================================== %
% ======================================================================= %













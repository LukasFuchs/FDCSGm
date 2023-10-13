clear
clc
% profile on
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\AdvectionProblem')
    addpath('..\..\StokesProblem')
    addpath('..\..\SetUp')
    addpath('..\..\ScaleParam')
else
    addpath('../../AdvectionProblem')
    addpath('../../StokesProblem')
    addpath('../../SetUp')
    addpath('../../ScaleParam')
end
% ======================================================================= %
T.tstart        =   tic;
%% ===================== Some initial definitions ======================= %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
Py.scale        =   'yes';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'none';
B.Aparam        =   'none';
B.DiffMethod    =   'none';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'variable';
B.EtaIni        =   'ellipse';
% inclusions bedingungen
B.RotAng        =   0;             % positive -> counter clockwise
B.EllA          =   2e2; % 1.75e3; [m]
B.EllB          =   2e2; % 0.25e3; [m]
B.T0            =   1000; 
B.TAmpl         =   1000; 
% ======================================================================= %
%% ================== Define initial temperature anomaly ================ %
B.Tini          =   'const';
Py.tparam       =   'const';
% ======================================================================= %
%% ========================= Define flow field ========================== %
B.IniFlow       =   'SimpleShear';
B.ebg           =   -1e-15;         % < 0 compression
B.FlowFac       =   [];
% ======================================================================= %
%% ==================== Define model geometry constants ================= %
M.H         =   -1;           %   Modeltiefe [ in km ]
M.xmax      =   1;               %   Seitenverhaeltniss
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   201;             %   Vertikale Gitteraufloesung
N.nx        =   201;             %   Horizontale Gitteraufloesung
% ======================================================================= %
%% =======================  Tracer advection ============================ %
N.nmx       =   5;
N.nmz       =   5;
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   10;                 %   Schwerebeschleunigung [m/s^2]
Py.rho0     =   3200;               %   Hintergunddichte [kg/m^3]
Py.k        =   3;                  %   Thermische Leitfaehigkeit [ W/m/K ]
Py.cp       =   1000;               %   Heat capacity [ J/kg/K ]
Py.alpha    =   5e-5;               %   Thermischer Expnasionskoef. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]

Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]

Py.eta0     =   1e23;           % Viskositaet [ Pa*s ]
Py.eta1     =   1e21;           % Inclusion viscosity
Py.rho1     =   Py.rho0;        % Includion density

Py.DeltaT   =   1000;           % Temperaturdifferenz
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Geschwindigkeitsrandbedingungen --------------------------------------- %
%   0 - no slip; 1 - free slip;
B.tbc       =   0;              %   Obenw
B.bbc       =   0;              %   Unten
B.lbc       =   0;              %   Links
B.rbc       =   0;              %   Rechts

% Thermische Randbedingungen -------------------------------------------- %
B.ttbc      =   'const';
B.btbc      =   'const';
B.ltbc      =   'const';
B.rtbc      =   'const';
% Define actual temperature condition: 
%   if 'flux'   -   following parameters define the flux through the boundary
%   if 'const'  -   following parameters define the temperatue at the boundary in K
% Bottom, e.g. 1e-3
B.lhf       =   1000;
B.rhf       =   1000;
B.thf       =   1000;
B.bhf       =   1000;
% ======================================================================= %
%% ====================== Define time constants ========================= %
T.tmaxini   =   4500;           %   Maximale Zeit in Ma
T.itmax     =   1;              %   Maximal erlaubte Anzahl der Iterationen
T.dtfac     =   1.0;            %   Advektionscourantkriterium
T.dtdifac   =   1.0;            %   Diffusions Stabilitaetskriterium
% ======================================================================= %
%% ========================= Define fields required ===================== %
[Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
% ======================================================================= %
%% ======================== Setup initial conditions ==================== %
[T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
% ======================================================================= %
%% ======================= Rayleigh number conditions =================== %
Py.Ra           =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.eta0/Py.kappa;
% ======================================================================= %
%% ========================= Plot parameter ============================= %
Pl.inc      =   min(N.nz/10,N.nx/10);
Pl.inc      =   round(Pl.inc);
% ======================================================================= %
%% ========================== Scale Parameters ========================== %
switch lower(Py.scale)
    case 'yes'
        [M,N,D,T,S]         =   ScaleParameters(B,M,Py,N,D,T);
end
% ======================================================================= %
switch B.IniFlow
    case 'SimpleShear'
        gr  =   1;
        er  =   0;
    case 'PureShear'
        gr  =   0;
        er  =   -1;
end
if B.EllA==B.EllB
    [ Vx_N,Vx_S,Vx_W,Vx_E,Vz_N,Vz_S,Vz_W,Vz_E, Pa, Vxa, Vza ] = ...
        Dani_Solution_vec(M.x-M.L/2,M.z-M.H/2,M.x1-M.L/2,M.z1-M.H/2,...
        (B.EllA)/(-M.H*1e3),Py.eta1/Py.eta0,N.nx1,N.nz1,...
        gr,er);
    D.Pa =    Pa'; D.Vxa =  Vxa'; D.Vza =   Vza';
    clear Vxa Vza Pa
else
    error('Inclusion must be circular!')
end
Vx_S(1)     = Vx_W(1);  Vx_S(end)     = Vx_E(1);
Vx_N(1)     = Vx_W(end);Vx_N(end)     = Vx_E(end);
D.vxi(1:end-1,1)    =   Vx_W; D.vxi(1:end-1,N.nx) =   Vx_E;
D.vxi(1,:)          =   Vx_S; D.vxi(N.nz1,:)      =   Vx_N;
D.vzi(1:end,1)      =   Vz_W; D.vzi(1:end,N.nx1)  =   Vz_E;
D.vzi(1,1:end-1)    =   Vz_S; D.vzi(N.nz,1:end-1) =   Vz_N;
% keyboard
% ======================================================================= %
%% ========================= Time loop ================================= %%
for it = 1:T.itmax
    switch Py.eparam
        case 'const'
            [D,A]       =   solveSECE_const_EtaSc(D,Py,N,B,A);
            if (it == 2)
                N.beenhere = 1;
            end
        case 'variable'
            [D,A]       =   solveSECESc(D,Py,N,B);
        otherwise
                error('Viscosity not difined! Check Py.eparam parameter.')
    end
    % =================================================================== %
    %% =========== Interpolate velocity onto the regular grid =========== %
    [ID]        =   InterpStaggered(D,ID,N,'velocity');
%     D.meanV(it) =   mean(ID.v,'all');   % Mittleregeschwindigkeit
    D.meanV(it) =   rms(ID.vx(:) + ID.vz(:));
    % =================================================================== %
    [ID]        =   GetStrainRate(ID,N);
    ID.tauII    =   ID.eII.*D.eta.*2;
    ID.psi      =   ID.eII.*ID.tauII;
    incind      =   log10(D.eta)==log10(Py.eta1);
    % =================================================================== %
    %% ========================== Plot data ============================= %
    Pl.time     =   '';
    Pl.xlab     =   '$$x$$';
    Pl.zlab     =   '$$z$$';
    
    if (mod(it,5)==0||it==1)
        figure(1) % ----------------------------------------------------- %
        clf
        ax1 = subplot(2,2,1);
        plotfield(log10(D.eta),M.X,M.Z,Pl,'pcolor',...
            '$$log_{10} (\ \eta\ )$$','quiver',ID.vx,ID.vz);
        colormap(ax1,flipud(Pl.lapaz))
        ax2 = subplot(2,2,2);
        plotfield(log10(ID.psi),M.X,M.Z,Pl,'pcolor',...
            '$$log_{10} (\ \psi\ )$$');
        colormap(ax2,Pl.imola)
        ax3 = subplot(2,2,3);
        plotfield(log10(ID.eII),M.X,M.Z,Pl,'pcolor',...
            '$$log_{10} (\ \dot\varepsilon_{II}\ )$$')
        colormap(ax3,Pl.batlowW)
        ax4 = subplot(2,2,4);
        plotfield(log10(ID.tauII),M.X,M.Z,Pl,'pcolor',...
            '$$log_{10} (\ \tau_{II}\ )$$')
        colormap(ax4,Pl.nuuk)
    end       
    % =================================================================== %
    
end
% Staggered grid coordinates
[M.xVx,M.zVx] = meshgrid(M.x,M.z1);
[M.xVz,M.zVz] = meshgrid(M.x1,M.z);

D.Pe        = abs(D.P(2:end,2:end)-D.Pa);
D.vxe       = abs(D.vx(1:end-1,:)-D.Vxa);
D.vze       = abs(D.vz(:,1:end-1)-D.Vza);

figure(2)
ax1 = subplot(3,3,1); plotfield(D.Vxa,M.xVx,M.zVx,Pl,...
    'pcolor','$$vx\ (analytical)$$')
colormap(ax1,Pl.imola)
ax2 = subplot(3,3,2); plotfield(D.vx(1:end-1,:),M.xVx,M.zVx,Pl,...
    'pcolor','$$vx\ (numerical)$$')
colormap(ax2,Pl.imola)
ax3 = subplot(3,3,3); plotfield(D.vxe,M.xVx,M.zVx,Pl,...
    'pcolor','$$vx\ (error)$$')
colormap(ax3,Pl.batlowW)
ax4 = subplot(3,3,4); plotfield(D.Vza,M.xVz,M.zVz,Pl,...
    'pcolor','$$vz\ (analytical)$$')
colormap(ax4,Pl.imola)
ax5 = subplot(3,3,5); plotfield(D.vz(:,1:end-1),M.xVz,M.zVz,Pl,...
    'pcolor','$$vz\ (numerical)$$')
colormap(ax5,Pl.imola)
ax6 = subplot(3,3,6); plotfield(D.vze,M.xVz,M.zVz,Pl,...
    'pcolor','$$vz\ (error)$$')
colormap(ax6,Pl.batlowW)
ax7 = subplot(3,3,7); plotfield(D.Pa,M.X1,M.Z1,Pl,...
    'pcolor','$$P\ (analytical)$$')
colormap(ax7,Pl.hawaii)
ax8 = subplot(3,3,8); plotfield(D.P(2:end,2:end),M.X1,M.Z1,Pl,...
    'pcolor','$$P\ (numerical)$$')
colormap(ax8,Pl.hawaii)
ax9 = subplot(3,3,9); plotfield(D.Pe,M.X1,M.Z1,Pl,...
    'pcolor','$$P\ (error)$$')
colormap(ax9,Pl.batlowW)

psiinc1 = mean(ID.psi(incind));
psiinc2 = -sum(ID.psi(incind))./sum(incind(:).*N.dx.*N.dz);
psiinc3 = sum(ID.psi(incind))./pi/B.EllA/B.EllB;

%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
%     rmpath('..\..\DiffusionProblem')
    rmpath('..\..\AdvectionProblem')
    rmpath('..\..\StokesProblem')
    rmpath('..\..\SetUp')
    rmpath('..\..\ScaleParam')
else
%     rmpath('../../DiffusionProblem')
    rmpath('../../AdvectionProblem')
    rmpath('../../StokesProblem')
    rmpath('../../SetUp')
    rmpath('../../ScaleParam')
end
% ======================================================================= %

% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %

function [ Vx_N,Vx_S,Vx_W,Vx_E,Vy_N,Vy_S,Vy_W,Vy_E,Pa,Vxa,Vya ] = ...
    Dani_Solution_vec(xv,yv,xc,yc,rad,mus_i,nx,ny,gr,er)
% -------------------------------------------------------------------------
% ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION
% BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
% Vectorised version by:
% Thibault Duretz, Ludovic Raess - Unil 2016
% -------------------------------------------------------------------------
% INPUT:
% gr  =  1;                       % Simple shear: gr=1, er=0
% er  =  0;                       % Strain rate
mm  =  1;                       % Viscosity of matrix
mc  = mus_i;
A   = mm.*(mc-mm)./(mc+mm);
i   = sqrt(-1);
%-------------------------
Pa  = zeros(nx  ,ny  );
Vxa = zeros(nx+1,ny  );
Vya = zeros(nx  ,ny+1);
% PRESSURE
[XC2, YC2]   = ndgrid(xc, yc);
Z            = XC2 + i*YC2;
PH           = zeros(nx,ny);
PH( XC2.^2 + YC2.^2 <= rad.^2 ) = 1;
P            = -2.*mm.*(mc-mm)./(mc+mm).*real(rad^2./Z.^2.*(i*gr+2*er));  % outside inclusion
Pa(PH==0)    = P(PH==0);
% Conforming Nodes --------------------------------------------------------
% VELOCITY X
[XV2, YC2]   = ndgrid(xv, yc);
Z            = XV2 + i*YC2;
PH           = zeros(nx+1,ny);
PH( XV2.^2 + YC2.^2 <= rad.^2 ) = 1;
V_tot        = (mm/(mc+mm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z; % inside inclusion
Vxa(PH==1)   = real(V_tot(PH==1));
phi_z        = -(i/2)*mm*gr*Z-(i*gr+2*er)*A*rad^2*Z.^(-1);  % outside inclusion
d_phi_z      = -(i/2)*mm*gr + (i*gr+2*er)*A*rad^2./Z.^2;
conj_d_phi_z = conj(d_phi_z);
psi_z        = (i*gr-2*er)*mm*Z-(i*gr+2*er)*A*rad^4*Z.^(-3);
conj_psi_z   = conj(psi_z);
V_tot        = (phi_z- Z.*conj_d_phi_z - conj_psi_z) / (2*mm);
Vxa(PH==0)   = real(V_tot(PH==0));
% VELOCITY Y
[XC2, YV2]   = ndgrid(xc, yv);
Z            = XC2 + i*YV2;
PH           = zeros(nx,ny+1);
PH( XC2.^2 + YV2.^2 <= rad.^2 ) = 1;
V_tot        = (mm/(mc+mm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z; % inside inclusion
Vya(PH==1)   = imag(V_tot(PH==1));
phi_z        = -(i/2)*mm*gr*Z-(i*gr+2*er)*A*rad^2*Z.^(-1);  % outside inclusion
d_phi_z      = -(i/2)*mm*gr + (i*gr+2*er)*A*rad^2./Z.^2;
conj_d_phi_z = conj(d_phi_z);
psi_z        = (i*gr-2*er)*mm*Z-(i*gr+2*er)*A*rad^4*Z.^(-3);
conj_psi_z   = conj(psi_z);
V_tot        = (phi_z- Z.*conj_d_phi_z - conj_psi_z) / (2*mm);
Vya(PH==0)   = imag(V_tot(PH==0));
% Get BC
Vx_W  = Vxa(1   , :  )';
Vx_E  = Vxa(end , :  )';
Vy_S  = Vya(:   ,1   );
Vy_N  = Vya(:   ,end );
% Non Conforming Nodes ----------------------------------------------------
Vx_NC = zeros(nx+1,ny+1);
Vy_NC = zeros(nx+1,ny+1);
Vx_S  = zeros(nx+1,1   );
Vx_N  = zeros(nx+1,1   );
Vy_W  = zeros(1   ,ny+1)';
Vy_E  = zeros(1   ,ny+1)';
% VELOCITY X & Y -
[XV2, YV2]   = ndgrid(xv, yv);
Z            = XV2 + i*YV2;
PH           = zeros(nx+1,ny+1);
PH( XV2.^2 + YV2.^2 <= rad.^2 ) = 1;
V_tot        = (mm/(mc+mm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z; % inside inclusion
Vx_NC(PH==1) = real(V_tot(PH==1));
Vy_NC(PH==1) = imag(V_tot(PH==1));
phi_z        = -(i/2)*mm*gr*Z-(i*gr+2*er)*A*rad^2*Z.^(-1);  % outside inclusion
d_phi_z      = -(i/2)*mm*gr + (i*gr+2*er)*A*rad^2./Z.^2;
conj_d_phi_z = conj(d_phi_z);
psi_z        = (i*gr-2*er)*mm*Z-(i*gr+2*er)*A*rad^4*Z.^(-3);
conj_psi_z   = conj(psi_z);
V_tot        = (phi_z- Z.*conj_d_phi_z - conj_psi_z) / (2*mm);
Vx_NC(PH==0) = real(V_tot(PH==0));
Vy_NC(PH==0) = imag(V_tot(PH==0));
% Get BC
Vx_S(2:end-1,1) = Vx_NC(2:end-1,1  );
Vx_N(2:end-1,1) = Vx_NC(2:end-1,end);
Vy_W(2:end-1,1) = Vy_NC(1  ,2:end-1)';
Vy_E(2:end-1,1) = Vy_NC(end,2:end-1)';
end

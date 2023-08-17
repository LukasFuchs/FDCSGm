clear
clc
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\DiffusionProblem')
    addpath('..\SetUp')
else
    addpath('../DiffusionProblem')
    addpath('../SetUp')
end
% ======================================================================= %
%% ===================== Some initial definitions ======================= %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
Py.scale        =   'no';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'none';
B.Aparam        =   'none';
B.DiffMethod    =   'none';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'none';
% ======================================================================= %
%% ================== Define initial temperature anomaly ================ %
B.Tini          =   'ContGeotherm';
Py.tparam       =   'variable';
% ======================================================================= %
%% ========================= Define flow field ========================== %
B.IniFlow       =   'none';
B.FlowFac       =   [];
% ======================================================================= %
%% ==================== Define model geometry constants ================= %
M.H         =   -200;           %   Modeltiefe [ in km ]
M.xmax      =   1;              %   Seitenverhaeltniss

M.zUC       =   -10e3; % -25e3;          % Depth of the upper crust
M.zLC       =   -35e3;          % Depth of the lower crust

T           =   [];
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   201;             %   Vertikale Gitteraufloesung
N.nx        =   201;             %   Horizontale Gitteraufloesung
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   9.81;               %   Schwerebeschleunigung [m/s^2]
Py.alpha    =   -999;               %   dummy
% Mantle constants -------------------------------------
Py.rho0     =   3000;               %   Hintergunddichte [kg/m^3]
Py.k0       =   2.3;                %   Thermische Leitfaehigkeit [ W/m/K ]
Py.Q0       =   2.3e-12*Py.rho0;    %   Heat productino rate per mass [ W/kg ]

% Constants for the upper crust ------------------------
Py.rhoUC    =   2700;               %   [ kg/m^3 ]
Py.kUC      =   3.0;                %   [ W/m/K ]
Py.QUC      =   617e-12*Py.rhoUC;   %   [ W/kg ]

% Constant for the lower crust -------------------------
Py.rhoLC    =   2900;               %   [ kg/m^3 ]
Py.kLC      =   2.0;                %   [ W/m/K ]
Py.QLC      =   43e-12*Py.rhoLC;    %   [ W/kg ]

% Bottom temperature -----------------------------------
Py.Tpot      =   1315;              %   Potential temperature [ C ]
Py.Ttop      =   0; 
Py.Tbot      =   (Py.Tpot+Py.Ttop) + 0.5*abs(M.H);
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Thermische Randbedingungen -------------------------------------------- %
B.ttbc      =   'const';
B.btbc      =   'flux';
B.ltbc      =   'flux';
B.rtbc      =   'flux';

% Waermerandbedinungen
% Falls 'flux' - definiert es den Fluss
% Falls 'const' (fuer top und bottom) - definiert es die Temperatur in K
% Bottom, e.g. 1e-3
B.lhf       =   0;
B.rhf       =   0;
B.thf       =   Py.Ttop + 273.15;
B.bhf       =   0.01; %Py.Tbot + 273.15;% 0.04;     [ W/m^2 ]
% ======================================================================= %
%% ----------------------- Felddefinitionen ----------------------------- %
% Falls eine Strukturvariable schon zuvor definiert wurde, muss diese hier
% der Funktion als Eingabeargument mitgegeben werden. Das Selbe gilt fuer
% aller weiteren Funktionen.
[Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
% ----------------------------------------------------------------------- %

%% ---------------- Definition der Anfangsbedingungen ------------------- %
[T,D,B,Ma,Py]   =   SetUpInitialConditions(T,D,Py,M,N,B);

M.UCind         =   M.Z>=M.zUC;
M.LCind         =   M.Z>=M.zLC&M.Z<M.zUC;
M.Mind          =   M.Z<=M.zLC;

D.k(:)          =   Py.k0;
D.k(M.UCind)    =   Py.kUC;
D.k(M.LCind)    =   Py.kLC;

D.rho(:)        =   Py.rho0;
D.rho(M.UCind)  =   Py.rhoUC;
D.rho(M.LCind)  =   Py.rhoLC;

D.Q(:)          =   Py.Q0;
D.Q(M.UCind)    =   Py.QUC;
D.Q(M.LCind)    =   Py.QLC;
% ----------------------------------------------------------------------- %
%% ============================ Plot Parameter ========================== %
Pl.inc      =   min(N.nz/10,N.nx/10);
Pl.inc      =   round(Pl.inc);
set(figure(1),'position',[55.4,125.8,1326.4,636.2]);
% ======================================================================= %
%% ============================== Diffusion ============================= %
%   Steady state -> Poisson equation!
[D.T]           =   SolvePoisson2Dvaryk(D.Q,N,B,D.k);
% ======================================================================= %
%% =========================== Plotting data ============================ %
Pl.xlab     =   'x [ km ]';
Pl.zlab     =   'z [ km ]';
Pl.time     =   [];

figure(1)
subplot(2,3,1)
plotfield(D.Q,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '\itQ \rm\bf')
subplot(2,3,2)
plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '\it\rho \rm\bf')
subplot(2,3,4)
plotfield(D.k,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '\itk \rm\bf')
subplot(2,3,5)
plotfield(D.T,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '\itT \rm\bf')
subplot(1,3,3)
plot(D.T(:,N.nx1/2)-273.15 ,M.z./1e3,'LineWidth',2)
xlabel('T [ C ]'); ylabel('Depth [km]')
set(gca,'FontWeight','Bold','LineWidth',2)

figure(2) 
subplot(3,1,1)
plot(D.T(:,N.nx1/2)-273.15 ,M.z./1e3,'LineWidth',2)
xlabel('T [ C ]'); ylabel('Depth [km]')
set(gca,'FontWeight','Bold','LineWidth',2)
subplot(3,1,2)
plot(D.Q(:,N.nx1/2)./D.rho(:,N.nx1/2),M.z./1e3,'LineWidth',2)
xlabel('H [ W/kg ]'); ylabel('Depth [km]')
set(gca,'FontWeight','Bold','LineWidth',2)
q           = zeros(N.nz-1,1);
for j=1:N.nz-1
    q(j) = -(D.k(j+1,N.nx1/2) + D.k(j,N.nx1/2))/2 * ...
        (D.T(j+1,N.nx1/2) - D.T(j,N.nx1/2))/N.dz; 
end
q(1,1)      = -D.k(1,N.nx1/2)*(D.T(2,N.nx1/2)-D.T(1,N.nx1/2))/N.dz; 
q(N.nz,1)   = -D.k(N.nz,N.nx1/2)*(D.T(N.nz,N.nx1/2)-D.T(N.nz-1,N.nx1/2))/N.dz; 
subplot(3,1,3)
plot(q.*1e3,M.z./1e3,'LineWidth',2)
xlabel('q [ mW/m^2 ]'); ylabel('Depth [km]')
set(gca,'FontWeight','Bold','LineWidth',2)
% ======================================================================= %
%% ============================= Save data ============================== %
% DATA    = [D.T(:,N.nx1/2)-273.15,(M.z'./1e3)];
% save('Continental_Geotherm_uc_lc_2D.txt','DATA','-ascii')
% ======================================================================= %
%% ======================= Clear path structure ========================= %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\DiffusionProblem')
    rmpath('..\SetUp')
else
    rmpath('../DiffusionProblem')
    rmpath('../SetUp')
end
% ======================================================================= %

% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


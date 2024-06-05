clc
clear

if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\DiffusionProblem\')
    addpath('..\Rheology\')    
else
    addpath('../DiffusionProblem/')
    addpath('../Rheology/')
end

%% Define model ========================================================= %
% Model Constants ------------------------------------------------------- %
M.H         =   -200e3;         % Model height [ m ]

M.Type      =   'oceanic';
M.zUC       =   -25e3;          % [ m ]
M.zLC       =   -35e3;          % [ m ]
% ======================================================================= %

%% Physical constants =================================================== %
Py.g        =   9.81;           % Gravitational acceleration [ m/s^2 ]
% ======================================================================= %

%% Rheology constants =================================================== %
D.eII       =   1e-16;          % Background strain rate [ 1/s ]
D.d         =   1e-3;           % Grain size [ m ]
R.nrheo     =   100;            % Number of iterations

% Effective Viscosity ---
R.pressure  =   'yes';          %
R.type      =   'HK03';         %
R.cetaeff   =   'harmonic';     % 
R.itfac     =   0.7;            % Strain Rate Damping Factor
% Plasticity ---
R.plast     =   'Drucker-Prager';
% Grain Size ---
R.Grains    =   'no';
G.GSE       =   'be09';         %
G.growth    =   'no';           % 
G.itfac     =   0.7;            % Grain Size Damping Factor
% ======================================================================= %

%% Rheological parameters =============================================== %
% Parameters from Schmalholz et al., 2009
% Upper Crust ----------------------------------------------------------- %
Py.rhoUC    =   2700;           % Density [ kg/m^3 ]
Py.cpUC     =   1050;           % Specific heat [ J/kg/K ]
Py.kUC      =   2.5;            % Thermal conductivity [ W/m/K ]
Py.HUC      =   1.4e-6/Py.rhoUC;% Volumetric radiogenic heat source [ W/kg ]
% Plasticity ---
Py.chUC     =   10e6;           % Cohesion [ Pa ]
Py.frUC     =   30;             % Friction angle

% Lower Crust ----------------------------------------------------------- %
Py.rhoLC    =   2900;           % Density [ kg/m^3 ]
Py.cpLC     =   1050;           % Specific heat [ J/kg/K ]
Py.kLC      =   2.1;            % Thermal conductivity [ W/m/K ]
Py.HLC      =   0.4e-6/Py.rhoLC;% Volumetric readiogenic heat source [ W/kg ]
% Plasticity ---
Py.chLC     =   10e6;           % Cohesion [ Pa ]
Py.frLC     =   30;             % Friction angle

% Mantle parameters ----------------------------------------------------- %
Py.rhoM     =   3300;           % Density [ kg/m^3 ]
Py.cpM      =   1050;           % Specific heat [ J/kg/K ]
Py.kM       =   3.0;            % Thermal conductivity [ W/m/K ]
Py.HM       =   0.0;            % Volumetric radiogenic heat source [ W/kg ]
% Plasticity ---
Py.chM      =   10e6;           % Cohesion [ Pa ]
Py.frM      =   30;             % Friction angle
% ======================================================================= %

%% Temperature parameters =============================================== %
T.K0        =   273.15;
T.Tpot      =   1350 + T.K0;    % Potential temperautre [ K ]
T.dTadi     =   0.0;            % Adiabatic temperature gradient [ K/km ]
T.T0        =   273.15;         % Surface temperature [ K ]
T.ubound    =   'const';        % 
T.utbf      =   0;              % c     =   -k/q
T.lbound    =   'const';        % 
T.ltbf      =   0;              % c     =   -k/q
% ======================================================================= %

%% Time parameters ====================================================== %
t.dtfac     =   1.0;                % Courant criterion
t.age       =   1500;                 % Age of the lithosphere [ Ma ]
t.tfac      =   (60*60*24*365.25);  % Seconds per year
t.age       =   t.age.*1e6*t.tfac;  % Age in seconds
% ======================================================================= %

%% Numerical constants ================================================== %
N.debug     =   0;
N.nz        =   501;                   % Number of grid points
N.dz        =   M.H/(N.nz - 1);         % Grid resolution
M.z         =   (0:N.dz:M.H)';          % Depth [ m ]
M.zc        =   (N.dz/2:N.dz:M.H-N.dz/2)';
M.UCind     =   M.z>=M.zUC;             % Upper crust index
M.LCind     =   M.z>=M.zLC&M.z<M.zUC;   % Lower crust index
M.Mind      =   M.z<=M.zLC;             % Mantle index
% ======================================================================= %

%% Global constants ===================================================== %
GC.RG       =   8.3144;                 % Gas constant [ J/mol/K ]
plotparam   =   1;
% ======================================================================= %

%% Setup Model ========================================================== %
Py.fr       =   zeros(N.nz,1);
Py.ch       =   zeros(N.nz,1);
D.tauy      =   zeros(N.nz,1);

switch lower(M.Type)
    case 'oceanic'
        Py.fr(:)    =   Py.frM;
        Py.ch(:)    =   Py.chM;
    case 'continental'
        Py.fr(M.UCind)  =   Py.frM;
        Py.ch(M.UCind)  =   Py.chM;
        Py.fr(M.LCind)  =   Py.frM;
        Py.ch(M.LCind)  =   Py.chM;
        Py.fr(M.Mind)   =   Py.frM;
        Py.ch(M.Mind)   =   Py.chM;
end
% ======================================================================= %

%% Get temperature profile ============================================== %
switch lower(M.Type)
    case 'oceanic'
        [T,Py]  =   OceanicGeotherm_1D(T,M,N,Py,t,plotparam);
    case 'continental'
        [T,Py]  =   ContinentalGeotherm_1D(T,M,N,Py,t,plotparam);
    otherwise
        error('Temperature profile not defined correctly!')
end
% ======================================================================= %

%% Lithostatic Pressure ================================================= %
D.Plith         =   zeros(N.nz,1);
if size(Py.rho,1)==1
        D.Plith     =   -Py.rho*Py.g*M.z;     % [ Pa ]
else
    for i = 2:N.nz
        D.Plith(i)  =   D.Plith(i-1) + ...
            Py.rho(i)*Py.g*(M.z(i-1)-M.z(i));  % [ Pa ]
    end
end
D.P     =   D.Plith;
% ======================================================================= %

%% Plasticity =========================================================== %
switch lower(R.plast)
    case 'drucker-prager'
        D.tauy  =   D.P.*sind(Py.fr) + cosd(Py.fr).*Py.ch;
        D.etay  =   D.tauy./2/D.eII;
end
% ======================================================================= %

%% Effective Viscosity ================================================== %
[R]     =   SetRheoParam(R,M,N);

switch lower(R.Grains)
    case 'steadystate'
        [G,R]   =   SetGSEParam(G,R); 
end

[D]     =   RHEOLOGY(D,G,R,T,GC,N,M);
% ======================================================================= %

%% Total strength of the lithosphere ==================================== %
% Integral of tauII over the entire depth profile ----------------------- %
D.S     =   -trapz(M.z,D.tauII);
% ======================================================================= %

%% Plot data ============================================================ %
figure
clf
subplot(1,4,1)
plot(T.T-T.K0,M.z./1e3,'r-','LineWidth',2)
xlabel('$$T [ ^{\circ}C ]$$','Interpreter','latex')
ylabel('$$Depth [ km ]$$','Interpreter','latex')
title([{'$$Temperature$$'},{['$$age: $$',num2str(t.age/1e6/t.tfac),'$$ Ma$$']}],...
    'Interpreter','latex')
set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,...
    'TickLabelInterpreter','latex')

subplot(1,4,2)
plot(D.eta,M.z./1e3,'k-','LineWidth',2)
hold on
plot(D.etaf,M.z./1e3,'b--','LineWidth',2)
plot(D.etal,M.z./1e3,'r--','LineWidth',2)
plot(D.etay,M.z./1e3,'m:','LineWidth',2)
legend('effectiv','diff','disl','yield','Location','SouthEast')
xlabel('$$\eta [ Pa s ]$$','Interpreter','latex')
ylabel('$$Depth [ km ]$$','Interpreter','latex')
title([{'$$Viscosity$$'},{['$$age: $$',num2str(t.age/1e6/t.tfac),'$$ Ma$$']}],...
    'Interpreter','latex')
axis([1e15 1e30 M.H./1e3 0])
set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,'xscale','log',...
    'TickLabelInterpreter','latex')

subplot(1,4,3)
plot(D.epstot,M.z./1e3,'k-','LineWidth',2)
hold on
plot(D.epsf,M.z./1e3,'b--','LineWidth',2)
plot(D.epsl,M.z./1e3,'r--','LineWidth',2)
legend('total','diff','disl','Location','SouthWest')
xlabel('$$\epsilon [ 1/s ]$$','Interpreter','latex')
ylabel('$$Depth [ km ]$$','Interpreter','latex')
title([{'$$Strain Rate$$'},{['$$age: $$',num2str(t.age/1e6/t.tfac),'$$ Ma$$']}],...
    'Interpreter','latex')
set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,'xscale','log',...
    'TickLabelInterpreter','latex')

subplot(1,4,4)
plot(D.tauII./1e6,M.z./1e3,'k-','LineWidth',2)
% hold on
% plot(D.tauy./1e6,M.z./1e3,'r--','LineWidth',2)
xlabel('$$\tau_{II} [ MPa ]$$','Interpreter','latex')
ylabel('$$Depth [ km ]$$','Interpreter','latex')
title([{'$$Stress$$'},{['$$age: $$',num2str(t.age/1e6/t.tfac),'$$ Ma$$']}],...
    'Interpreter','latex')
set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,'xscale','lin',...
    'TickLabelInterpreter','latex')

if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\DiffusionProblem\')
    rmpath('..\Rheology\')    
else
    rmpath('../DiffusionProblem/')
    rmpath('../Rheology/')
end

% ======================================================================= %
% ================================ END ================================== %
% ======================================================================= %
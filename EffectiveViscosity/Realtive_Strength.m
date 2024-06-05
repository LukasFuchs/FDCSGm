clc
clear

if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\DiffusionProblem\')
    addpath('..\Rheology\')
    addpath('..\SetUp\')
else
    addpath('../DiffusionProblem/')
    addpath('../Rheology/')
    addpath('../SetUp/')
end

Type    =   {'oceanic','continental'};
ebg     =   1e-15;
nn      =   30;
tminO   =   30;
tmaxO   =   200;
ageO    =   linspace(tminO,tmaxO,nn);
tminC   =   500;
tmaxC   =   2000;
ageC    =   linspace(tminC,tmaxC,nn);
age     =   [ageO;ageC];
S       =   zeros(size(Type,2),nn);
plotparam   =   1;

Pl      =   [];

[~,~,~,~,~,~,~,Pl]  =   SetUpFields([],[],[],[],[],Pl);

%% Define model ========================================================= %
% Model Constants ------------------------------------------------------- %
M.H         =   -410e3;     % Model height [ m ]
M.zUC       =   -25e3;      % [ m ]
M.zLC       =   -35e3;      % [ m ]
% ======================================================================= %
%% Physical constants =================================================== %
Py.g        =   9.81;       % Gravitational acceleration [ m/s^2 ]
% ======================================================================= %
%% Rheology constants =================================================== %
D.eII       =   ebg;        % Background strain rate [ 1/s ]
D.d         =   1e-3;       % Grain size [ m ]
R.nrheo     =   100;        % Number of iterations
% Effective Viscosity ---
R.pressure  =   'yes';      %
R.type      =   'HK03';     %
R.cetaeff   =   'harmonic'; %
R.itfac     =   0.7;        % Strain Rate Damping Factor
% Plasticity ---
R.plast     =   'drucker-prager';      %
% Grain Size ---
R.Grains    =   'no';
G.GSE       =   'be09';     %
G.growth    =   'no';       %
G.itfac     =   0.7;        % Grain Size Damping Factor
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
T.Tpot      =   1350 + T.K0;% Potential temperautre [ K ]
T.dTadi     =   0.0;        % Adiabatic temperature gradient [ K/km ]
T.T0        =   273.15;     % Surface temperature [ K ]
T.ubound    =   'const';    %
T.utbf      =   0;          % c     =   -k/q
T.lbound    =   'const';    %
T.ltbf      =   0;          % c     =   -k/q
% ======================================================================= %
%% Time parameters ====================================================== %
t.dtfac     =   1.0;        % Courant criterion
t.tfac      =   (60*60*24*365.25);  % Seconds per year
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
GC.RG       =   8.3144;         % Gas constant [ J/mol/K ]
% ======================================================================= %
%% LOOPS ================================================================ %
for k = 1:size(Type,2)
    disp(Type{k})
    %% Define model ===================================================== %
    % Model Constants --------------------------------------------------- %
    M.Type      =   Type{k};
    % =================================================================== %
    %% Setup Model ====================================================== %
    Py.fr       =   zeros(N.nz,1);
    Py.ch       =   zeros(N.nz,1);
    D.tauy      =   zeros(N.nz,1);
    switch lower(M.Type)
        case 'oceanic'
            Py.fr(:)        =   Py.frM;
            Py.ch(:)        =   Py.chM;
        case 'continental'
            Py.fr(M.UCind)  =   Py.frM;
            Py.ch(M.UCind)  =   Py.chM;
            Py.fr(M.LCind)  =   Py.frM;
            Py.ch(M.LCind)  =   Py.chM;
            Py.fr(M.Mind)   =   Py.frM;
            Py.ch(M.Mind)   =   Py.chM;
    end
    % =================================================================== %
    %% Effective Viscosity ============================================== %
    [R]     =   SetRheoParam(R,M,N);
    switch lower(R.Grains)
        case 'steadystate'
            [G,R]   =   SetGSEParam(G,R);
    end
    for l = 1:nn
        disp(age(k,l))
        %% Time parameters ============================================== %
        t.age       =   age(k,l);   % Age of the lithosphere [ Ma ]
        t.age       =   t.age.*1e6*t.tfac;  % Age in seconds
        % =============================================================== %
        %% Get temperature profile ====================================== %
        switch lower(M.Type)
            case 'oceanic'
                [T,Py]  =   OceanicGeotherm_1D(T,M,N,Py,t,0);
            case 'continental'
                [T,Py]  =   ContinentalGeotherm_1D(T,M,N,Py,t,0);
            otherwise
                error('Temperature profile not defined correctly!')
        end
        % =============================================================== %
        %% Lithostatic Pressure ========================================= %
        D.Plith         =   zeros(N.nz,1);
        if size(Py.rho,1)==1
            D.Plith     =   -Py.rho*Py.g*M.z;           % [ Pa ]
        else
            for i = 2:N.nz
                D.Plith(i)  =   D.Plith(i-1) + ...
                    Py.rho(i)*Py.g*(M.z(i-1)-M.z(i));   % [ Pa ]
            end
        end
        D.P     =   D.Plith;
        % =============================================================== %
        %% Plasticity =================================================== %
        switch lower(R.plast)
            case 'drucker-prager'
                D.tauy  =   D.P.*sind(Py.fr) + cosd(Py.fr).*Py.ch;
                D.etay  =   D.tauy./2/D.eII;
        end
        % =============================================================== %
        %% Effective Viscosity ========================================== %
        [D]     =   RHEOLOGY(D,G,R,T,GC,N,M);
        % =============================================================== %
        %% Total strength of the lithosphere ============================ %
        % Integral of tauII over the entire depth profile --------------- %
        D.S     =   -trapz(M.z,D.tauII);
        % =============================================================== %
        S(k,l)  =   D.S;
        %% Plot data ==================================================== %
        if plotparam==1
            set(figure(2),'Position',[1.8,1.8,766.4,780.8])
            figure(2)
            clf
            subplot(1,2,1)
            plot(D.tauII./1e6,M.z./1e3,'k-','LineWidth',2)
            hold on
            plot(T.T,M.z./1e3,'r--','LineWidth',2);
            text(100,-270,'$$\tau_{II}\ [ MPa ]\ /$$','Interpreter','latex',...
                'FontSize',20,'FontWeight','bold')
            text(1450,-270,'$$T\ [ K ]$$','Interpreter','latex',...
                'FontSize',20,'FontWeight','bold','Color','r')
            ylabel('$$Depth [ km ]$$','Interpreter','latex')
            axis([0 2000 M.H/1e3 0])
            text(100,-230,[{['$$ t =\ $$',...
                sprintf('%1.2f',(t.age/1e6/t.tfac)),'$$ Ma$$']},...
                {['$$ \dot{\varepsilon}_{bg} =\ $$',sprintf('%4.0e',D.eII),...
                '$$\ s^{-1}$$']}],...
                'Interpreter','latex','FontSize',15,'FontWeight','bold')
            set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',20,...
                'xscale','lin',...
                'TickLabelInterpreter','latex')
            
            subplot(1,2,2)
            plot(D.eta,M.z./1e3,'k-','LineWidth',2)
            hold on
            plot(D.etaf,M.z./1e3,'b--','LineWidth',2)
            plot(D.etal,M.z./1e3,'r--','LineWidth',2)
            plot(D.etay,M.z./1e3,'m:','LineWidth',2)
            legend('effectiv','diff','disl','yield','Location','SouthEast')
            xlabel('$$\eta [ Pa s ]$$','Interpreter','latex')
            ylabel('$$Depth [ km ]$$','Interpreter','latex')
            %             title([{'$$Viscosity$$'},{['$$age: $$',num2str(t.age/1e6/t.tfac),'$$ Ma$$']}],...
            %                 'Interpreter','latex')
            axis([1e17 1e26 M.H./1e3 0])
            set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',20,...
                'xscale','log',...
                'TickLabelInterpreter','latex')
            
            saveas(figure(2),['data/Plot_',M.Type,'_',...
                num2str(round(age(k,l)))],'png')
            saveas(figure(2),['data/Plot_',M.Type,'_',...
                num2str(round(age(k,l)))],'svg')
            saveas(figure(2),['data/Plot_',M.Type,'_',...
                num2str(round(age(k,l)))],'fig')
        end
    end
end
% Relative strength between oceanic and continental lithosphere --------- %
% RS    =   ( (S_CC - S_OC) / S_CC ) * 100
% RS      =   (1./S(2,:)')*(S(2,:)-S(1,:)).*100;
mark    =   {'w-','w-','r:','r:','r:','y--','y--','y--','y--'};
fall    =   1;
wo      =   'UM';
switch fall
    case 1
        load(['data\CLoOL_',wo,'.mat']);
    case 2
        load(['data\DS_CL_',wo,'.mat']);
    case 3
        load(['data\DS_OL_',wo,'.mat']);
end
RS  = zeros(size(Type,2),nn);
for k = 1:nn
    for l = 1:nn
        switch fall
            case 1
                RS(l,k)     =   ( S(2,l) / S(1,k) );
            case 2
                RS(l,k)     =   ( (S(2,l)-S(1,k)) / S(2,l) )*100;
            case 3
                RS(l,k)     =   ( (S(2,l)-S(1,k)) / S(1,k) )*100;
        end
    end
end
[XX,YY] =   meshgrid(ageO,ageC);
set(figure(1),'Position',[279.4,107.4,869.6,584])
ax = figure(1);
clf
contourf(XX,YY,RS,'w--','ShowText','on');
hold on
for j = 1:length(Data(:,1))
    contour(XX,YY,RS,[Data(j,1) Data(j,1)],[mark{j}],'ShowText','on');
end
colormap(ax,Pl.batlowW)
xlabel('$$t_{Oceanic}\ [Ma]$$','Interpreter','latex');
ylabel('$$t_{Continental}\ [Ma]$$','Interpreter','latex')
title('$$Relative\ Strength\ S_R\ =\ \frac{S_{CL}}{S_{OL}},\ where\ S_i = \int_{0}^{-H} \tau_{II,i}\ dz$$','Interpreter','latex')
cb = colorbar; % shading interp, lighting phong
switch fall
    case 1
        title(cb,'factor','Interpreter','latex')
    case 2
        title(cb,'%','Interpreter','latex')
    case 3
        title(cb,'%','Interpreter','latex')
end
axis normal; box on
set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,'xscale','log',...
    'yscale','log','TickLabelInterpreter','latex')

DATA = {XX YY age RS S};
save(['data/RelativeStrength_',num2str(fall),'_',wo,'.mat'],'DATA','-mat')
saveas(figure(1),['data/RelativeStrength_',num2str(fall),'_',wo],'png')
saveas(figure(1),['data/RelativeStrength_',num2str(fall),'_',wo],'svg')
saveas(figure(1),['data/RelativeStrength_',num2str(fall),'_',wo],'fig')

if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\DiffusionProblem\')
    rmpath('..\Rheology\')
    rmpath('..\SetUp\')
else
    rmpath('../DiffusionProblem/')
    rmpath('../Rheology/')
    rmpath('../SetUp/')
end
% ======================================================================= %
% ================================ END ================================== %
% ======================================================================= %
clear
clc
% profile on
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\DiffusionProblem')
    addpath('..\..\AdvectionProblem')
    addpath('..\..\StokesProblem')
    addpath('..\..\SetUp')
    addpath('..\..\ScaleParam')
    % Load benchmark data
    Ger     =   load('data\Gerya2019.txt');
else
    addpath('../../DiffusionProblem')
    addpath('../../AdvectionProblem')
    addpath('../../StokesProblem')
    addpath('../../SetUp')
    addpath('../../ScaleParam')
    % Load benchmark data
    Ger     =   load('data/Gerya2019.txt');
end
% ======================================================================= %
%% Some initial definitions --------------------------------------------- %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
Py.scale        =   'yes';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'semi-lag';
B.Aparam        =   'temp';
B.DiffMethod    =   'explicit';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'const';
Py.b            =   log(1000); %log(16384);          % Temperaturabhaengigkeit
Py.c            =   0; %log(64);                  % Tiefenabhaengigkeit
B.EtaIni        =   'tdep';
% ======================================================================= %
%% ================== Define initial temperature anomaly ================ %
B.Tini          =   'block';
B.T0            =   1000;
B.TAmpl         =   20; 
Py.tparam       =   'const';
% ======================================================================= %
%% ========================= Define flow field ========================== %
% B.IniFlow       =   'none';
B.FlowFac       =   10;
% ======================================================================= %
%% ==================== Define model geometry constants ================= %
M.H         =   -1000;          %   Modeltiefe [ in km ]
M.xmax      =   1;              %   Seitenverhaeltniss
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   21;             %   Vertikale Gitteraufloesung
N.nx        =   21;             %   Horizontale Gitteraufloesung
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   10;                 %   Schwerebeschleunigung [m/s^2]
Py.rho0     =   4000;               %   Hintergunddichte [kg/m^3]
Py.k        =   5;                  %   Thermische Leitfaehigkeit [ W/m/K ]
Py.cp       =   1250;               %   Heat capacity [ J/kg/K ]
Py.alpha    =   2.5e-5;             %   Thermischer Expnasionskoef. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]

Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]

Py.eta0     =   1e23;           %   Viskositaet [ Pa*s ]

Py.DeltaT   =   1000;           % Temperaturdifferenz

% Falls Ra < 0 gesetzt ist, dann wird Ra aus den obigen Parametern
% berechnet. Falls Ra gegeben ist, dann wird die Referenzviskositaet so
% angepasst, dass die Skalierungsparameter die gegebene Rayleigh-Zahl
% ergeben.
Py.Ra       =   -999;           % Rayleigh number
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Geschwindigkeitsrandbedingungen --------------------------------------- %
%   0 - no slip; 1 - free slip;
B.tbc       =   1;              %   Oben
B.bbc       =   1;              %   Unten
B.lbc       =   1;              %   Links
B.rbc       =   1;              %   Rechts

% Thermische Randbedingungen -------------------------------------------- %
B.ttbc      =   'const';
B.btbc      =   'const';
B.ltbc      =   'flux';
B.rtbc      =   'flux';
% Define actual temperature condition: 
%   if 'flux'   -   following parameters define the flux through the boundary
%   if 'const'  -   following parameters define the temperatue at the boundary in K
B.lhf       =   0;
B.rhf       =   0;
B.thf       =   273;
B.bhf       =   B.thf + Py.DeltaT;
% ======================================================================= %
%% ====================== Define time constants ========================= %
T.tmaxini   =   10000;          %   Maximale Zeit in Ma
T.itmax     =   1e6;            %   Maximal erlaubte Anzahl der Iterationen
T.dtfac     =   1.0;            %   Advektionscourantkriterium
T.dtdifac   =   0.9;            %   Diffusions Stabilitaetskriterium
% ======================================================================= %
%% ======================= Rayleigh number conditions =================== %
if Py.Ra < 0
    % Falls die Rayleigh Zahl nicht explizit angegeben wird, wird sie hier
    % berechnet.
    Py.Ra   =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H*1e3)^3/Py.eta0/Py.kappa;
else
    % Falls die Rayleigh Zahl gegeben ist müssen wir eine Variable
    % anpassen, z.B. die Referenzviskosität.
    Py.eta0 =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H*1e3)^3/Py.Ra/Py.kappa;
end
% ======================================================================= %
% n           =   [2,3,4,5];
n           =   [2 3];
nz          =   ceil((n-1).*(Py.Ra/4)^(1/3)+1);
Nus         =   zeros(length(n),1);
VRMS        =   zeros(length(n),1);
for i = 1:length(n)
    T.tstart    =   tic;
    N.nz    =   nz(i);
    N.nx    =   nz(i);
    % =================================================================== %
    %% ===================== Define fields required ===================== %
    [Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
    % =================================================================== %
    %% ==================== Setup initial conditions ==================== %
    [T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
    % =================================================================== %
    %% ===================== Plot parameter ============================= %
    Pl.inc      =   min(N.nz/10,N.nx/10);
    Pl.inc      =   round(Pl.inc);
    Pl.xlab     =   '$$x$$';
    Pl.zlab     =   '$$z$$';
    switch Pl.plotfields
        case 'yes'
            if strcmp(getenv('OS'),'Windows_NT')
                set(figure(1),'position',[1.8,1.8,766.4,780.8]);
                Pl.h        =   figure(1);
            else
                set(figure(1),'position',[-1919,1,960,988]);
                Pl.h        =   figure(1);
            end
    end
    % Animation settings ------------------------------------------------ %
    switch Pl.savefig
        case 'yes'
            M.ModDir    = ['data/resolution_test/Blanckenbach_Ra_',...
                sprintf('%2.2e',Py.Ra),...
                '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
                '_nx_',num2str(N.nz),'_nz_',num2str(N.nz)];
            if ~exist(M.ModDir,'dir')
                mkdir(M.ModDir)
            end
            Pl.filename     =   [M.ModDir,'/Evolution.gif'];
            set(figure(1),'position',[1.8,1.8,766.4,780.8]);
            Pl.h            =   figure(1);
    end
    % =================================================================== %
    %% ====================== Scale Parameters ========================== %
    switch lower(Py.scale)
        case 'yes'
            [M,N,D,T,S]         =   ScaleParameters(B,M,Py,N,D,T);
    end
    % =================================================================== %
    %% ============ Information for the command window ================== %
    fprintf([' Thermische Konvektion fuer Ra = %2.2e:\n  --------------------- ',...
        '\n Diffusion mit: %s',...
        '\n Advektion mit: %s',...
        '\n Viskositaet ist: %s',...
        '\n Anfangstemperaturfeld: %s',...
        '\n Aufloesung (nx x nz): %i x %i',...
        '\n Refrenzviskositaet [Pa s]: %2.2e',...
        '\n  --------------------- ',...
        '\n '],Py.Ra,B.DiffMethod,B.AdvMethod,Py.eparam,B.Tini,...
        N.nx,N.nz,Py.eta0);
    fprintf(['Maximum Time : %1.4g',...
        '\n Maximale Anzahl an Iterationen: %5i',...
        '\n  --------------------- \n'],...
        T.tmax,T.itmax)
    % =================================================================== %
    %% ===================== Time loop ================================= %%
    for it = 1:T.itmax
        if(strcmp(B.AdvMethod,'none')==0)
            switch Py.eparam
                case 'const'
                    [D,A]       =   solveSECE_const_EtaSc(D,Py,N,B,A);
                    if (it == 2)
                        N.beenhere = 1;
                    end
                case 'variable'
                    [D,A]       =   solveSECESc(D,Py,N,B);
                otherwise
                    error('Viskositaet nicht definiert! Siehe Parameter Py.eparam.')
            end
        end
        % =============================================================== %
        %% ================= Calculate time stepping ==================== %
        T.dt        =   T.dtfac*min(N.dx,abs(N.dz))/...
            (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));
        T.dtdiff    =   T.dtdifac*min([N.dx,abs(N.dz)])^2/4;
        T.dt        =   min(T.dt,T.dtdiff);
        if it>1
            T.time(it)  =   T.time(it-1) + T.dt;
            if T.time(it) > T.tmax
                T.dt        =   T.tmax - T.time(it-1); 
                T.time(it)  =   T.time(it-1) + T.dt;
            end
        end
        % =============================================================== %
        %% ======= Interpolate velocity onto the regular grid =========== %
        [ID]        =   InterpStaggered(D,ID,N,'velocity');
        %     D.meanV(it) = mean(ID.v,'all');   % Mittleregeschwindigkeit
        D.meanV(it) = rms(ID.vx(:) + ID.vz(:));
        % =============================================================== %
        %% ====================== Plot data ============================= %
        Pl              =   PlotData(it,Pl,T,D,M,ID,Py);
        % =============================================================== %
        %% ====================== Advection ============================= %
        [D,Ma,ID]       =   Advection(it,N,B,D,ID,Py,T.dt,M,Ma);
        % =============================================================== %
        %% ====================== Diffusion ============================= %
        D.T             =   Diffusion(B,D.T,D.Q,D.rho,Py,T.dt,N);
        % =============================================================== %
        %% ============== Heat flow at the surface ====================== %
        D.dTtop         =   (D.T(1,:)-D.T(2,:))./N.dz;
        D.dTbot         =   (D.T(end-1,:)-D.T(end,:))./N.dz;
        D.Nus(it)       =   mean(D.dTtop);

        D.meanT(:,it)   =   mean(D.T,2);
        % =============================================================== %
        %% ================== Update viscosity ========================== %
        switch lower(B.EtaIni)
            case 'tdep'
                D.eta   =   exp( -Py.b.*((D.T-D.T(1,1))./(D.T(N.nz,1)-D.T(1,1)))...
                    + Py.c.*M.Z);
        end
        % =============================================================== %
        %% ====================== Check break =========================== %
        [answer,T]  =   CheckBreakCriteria(it,T,D,M,Pl,ID,Py);
        switch answer
            case 'yes'
                break
        end
    end
    switch Pl.savefig
        case 'yes'
            saveas(figure(1),[M.ModDir,'/Field_SS'],'png')
    end
    % =================================================================== %
    %% ==================== Plot time serieses ========================== %
    PlotTimeSerieses(Py,T,D,M,N,Ger)
    switch Pl.savefig
        case 'yes'
            saveas(figure(2),[M.ModDir,'/TimeSeries'],'png')
            save([M.ModDir,'/B.mat'],'B','-mat')
            save([M.ModDir,'/D.mat'],'D','-mat')
            save([M.ModDir,'/M.mat'],'M','-mat')
            save([M.ModDir,'/N.mat'],'N','-mat')
            save([M.ModDir,'/Pl.mat'],'Pl','-mat')
            save([M.ModDir,'/Py.mat'],'Py','-mat')
            save([M.ModDir,'/S.mat'],'S','-mat')
            save([M.ModDir,'/T.mat'],'T','-mat')
    end
    % =================================================================== %
    T.tend      =   toc(T.tstart);

    Nus(i)      =   mean(D.Nus(T.tind:T.indtime));
    VRMS(i)     =   mean(D.meanV(T.tind:T.indtime));

end
% ======================================================================= %
if (Py.eta0 == 1e23)
    blmod   = 1;
elseif (Py.eta0 == 1e22)
    blmod   = 2;
elseif (Py.eta0 == 1e21)
    blmod   = 3;
end

figure(3)
subplot(2,1,1)
plot(1./nz.^2,Nus,'o','MarkerSize',8,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r')
hold on
plot(1e-6,Ger(1,blmod),'s','MarkerSize',8,'MarkerEdgeColor','r',...
    'MarkerFaceColor','k')
set(gca,'FontWeight','Bold','yscale','log','xscale','log',...
    'LineWidth',2,'FontSize',15,'TickLabelInterpreter','latex')
xlabel('$$1/nz/nx$$','Interpreter','latex')
ylabel('$$Nusselt\ Number$$','Interpreter','latex')
title('$$Resolution\ test$$','Interpreter','latex')

subplot(2,1,2)
plot(1./nz.^2,VRMS,'o','MarkerSize',8,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r')
hold on
plot(1e-6,Ger(2,blmod),'s','MarkerSize',8,'MarkerEdgeColor','r',...
    'MarkerFaceColor','k')
set(gca,'FontWeight','Bold','yscale','log','xscale','log',...
    'LineWidth',2,'FontSize',15,'TickLabelInterpreter','latex')
xlabel('$$1/nz/nx$$','Interpreter','latex')
ylabel('$$V_{RMS}$$','Interpreter','latex')
legend('$$Model$$','$$Benchmark$$','Interpreter','latex')

DATA    =   [nz',Nus,VRMS];
switch Pl.savefig
    case 'yes'
        if ~exist('data/resolution_test/','dir')
            mkdir('data/resolution_test/')
        end
        save(['data/resolution_test/ResTestData',num2str(Py.eta0),'.mat'],'DATA')
        saveas(figure(3),['data/resolution_test/ResTest_eta_',Py.eparam,...
            ],'png')
        saveas(figure(3),['data/resolution_test/ResTest_eta_',Py.eparam,...
            ],'fig')
end

%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\..\DiffusionProblem')
    rmpath('..\..\AdvectionProblem')
    rmpath('..\..\StokesProblem')
    rmpath('..\..\SetUp')
    rmpath('..\..\ScaleParam')
else
    rmpath('../../DiffusionProblem')
    rmpath('../../AdvectionProblem')
    rmpath('../../StokesProblem')
    rmpath('../../SetUp')
    rmpath('../../ScaleParam')
end
% ======================================================================= %

% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


clear
clc

%% Add some paths ------------------------------------------------------- %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\DiffusionProblem')
    addpath('..\AdvectionProblem')
    addpath('..\StokesProblem')
    addpath('..\SetUp')
else
    addpath('../DiffusionProblem')
    addpath('../AdvectionProblem')
    addpath('../StokesProblem')
    addpath('../SetUp')
end

T.tstart        =   tic;

%% Some initial definitions --------------------------------------------- %
% Define method for solving the energy equation ------------------------- %
B.AdvMethod     =   'semi-lag';
B.DiffMethod    =   'ADI';

% Variable oder konstante Vikositaet ------------------------------------ %
Py.eparam       =   'const';
Py.b            =   log(1000);          % Temperaturabhaengigkeit
Py.c            =   0;                  % Tiefenabhaengigkeit

% Define initial temperature anomaly ------------------------------------ %
B.Tini          =   'ellipse';

% Define flow field ----------------------------------------------------- %
B.IniFlow       =   'none';
B.FlowFac       =   10;
% ----------------------------------------------------------------------- %

%% --------------------- Definition der Modelkonstanten ----------------- %
M.H         =   -2971;          %   Modeltiefe [ in km ]
M.xmax      =   3;              %   Seitenverhaeltniss
% ----------------------------------------------------------------------- %

%% ------------------- Definition des Numerischen Gitters --------------- %
N.nz        =   51;             %   Vertikale Gitteraufloesung
N.nx        =   151;             %   Horizontale Gitteraufloesung
% ----------------------------------------------------------------------- %

%% -------------- Definition Physikalischer Konstanten ------------------ %
Py.g        =   10;                 %   Schwerebeschleunigung [m/s^2]
Py.rho0     =   4000;               %   Hintergunddichte [kg/m^3]
Py.k        =   5;                  %   Thermische Leitfaehigkeit [ W/m/K ]
Py.cp       =   1250;               %   Heat capacity [ J/kg/K ]
Py.alpha    =   2.5e-5;             %   Thermischer Expnasionskoef. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]

Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]

Py.eta0     =   1e23;           %   Viskositaet [ Pa*s ]

Py.DeltaT   =   2500;           % Temperaturdifferenz

Py.Ra       =   1e6;           % Rayleigh number
% ----------------------------------------------------------------------- %

%% -------------------- Definition der Randbedingungen ------------------ %
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

% Waermerandbedinungen
% Falls 'flux' - definiert es den Fluss
% Falls 'const' (fuer top und bottom) - definiert es die Temperatur in K
% Bottom, e.g. 1e-3
B.lhf       =   0;
B.rhf       =   0;
B.thf       =   273;
B.bhf       =   B.thf + Py.DeltaT;
% ----------------------------------------------------------------------- %

%% ------------------------ Definition der Zeitkonstanten --------------- %
T.tmaxini   =   10000;          %   Maximale Zeit in Ma
T.itmax     =   1000;           %   Maximal erlaubte Anzahl der Iterationen
T.dtfac     =   10.0;           %   Advektionscourantkriterium
T.dtdifac   =   1.0;            %   Diffusions Stabilitaetskriterium

% T.time      =   zeros(T.itmax,1);
% T.tmax      =   T.tmax.*1e6*(365.25*24*60*60);      % in [ s ]

% N.beenhere  =   0;
% N.beenhere2 =   0;
% ----------------------------------------------------------------------- %

%% ----------------------- Felddefinitionen ----------------------------- %
% Falls eine Strukturvariable schon zuvor definiert wurde, muss diese hier
% der Funktion als Eingabeargument mitgegeben werden. Das Selbe gilt fuer
% aller weiteren Funktionen.
[Py,D,ID,M,N,A]     =   SetUpFields(Py,N,M,T);
% ----------------------------------------------------------------------- %

%% ---------------- Definition der Anfangsbedingungen ------------------- %
[D,B,Ma,Py]     =   SetUpInitialConditions(D,Py,M,N,B);

if Py.Ra < 0
    % Falls die Rayleigh Zahl nicht explizit angegeben wird, wird sie hier
    % berechnet.
    Py.Ra   =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.eta0/Py.kappa;
else
    % Falls die Rayleigh Zahl gegeben ist müssen wir eine Variable
    % anpassen, z.B. die Referenzviskosität.
    Py.eta0 =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.Ra/Py.kappa;
end
% ----------------------------------------------------------------------- %

%% ------------------------- Plot Parameter ----------------------------- %
Pl.inc      =   min(N.nz/10,N.nx/10);
Pl.inc      =   round(Pl.inc);

% Animation settings ---------------------------------------------------- %
filename    = ['Thermal_Conv_Ra_',sprintf('%2.2e',Py.Ra),...
    '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),'.gif'];
set(figure(1),'position',[1.8,1.8,766.4,780.8]);
h           =   figure(1);
% ----------------------------------------------------------------------- %

%% Scale Parameters ===================================================== %
[M,N,D,T,S]  = ScaleParameters(M,Py,N,D,T);
% ----------------------------------------------------------------------- %

% Information fuer das Befehlsfenster ----------------------------------- %
fprintf([' Thermische Konvektion fuer Ra = %2.2e:\n  --------------------- ',...
    '\n Diffusion mit: %s',...
    '\n Advektion mit: %s',...
    '\n Viskositaet ist: %s',...
    '\n Anfangstemperaturfeld: %s',...
    '\n Anfangsgeschwindigkeitsfeld: %s',...
    '\n  --------------------- ',...
    '\n '],Py.Ra,B.DiffMethod,B.AdvMethod,Py.eparam,B.Tini,B.IniFlow);
fprintf(['Maximum Time : %1.4g',...
    '\n Maximale Anzahl an Iterationen: %5i',...
    '\n  --------------------- \n'],...
    T.tmax,T.itmax)
% pause
% ----------------------------------------------------------------------- %

%% BEGINN DER ZEITSCHLEIFE ============================================= %%
for it = 1:T.itmax
    %% Erstellung der Koeffizienten Matrix und rechten Seite des
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
    % ------------------------------------------------------------------- %
    
    %% Berechnung der Zeitschrittlaenge --------------------------------- %
    T.dt        =   T.dtfac*min(N.dx,abs(N.dz))/...
        (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));
    
    T.dtdiff    =   T.dtdifac*min([N.dx,abs(N.dz)])^2/4;
    
    T.dt        =   min(T.dt,T.dtdiff);
    
    if it>1
        T.time(it)  =   T.time(it-1) + T.dt;
    end
    % ------------------------------------------------------------------- %
    
    %% Interpolation der Geschwindigkeiten auf das regulaere Gitter ----- %
    [ID]        =   InterpStaggered(D,ID,N,'velocity');
    D.meanV(it) = mean(ID.v,'all');   % Mittleregeschwindigkeit
    % ------------------------------------------------------------------- %
    
    %% Darstellung der Daten -------------------------------------------- %
    Pl.time     =   ...
        ['@ Iteration: ',sprintf('%i',it),...
        '; Time: ',sprintf('%2.2e',T.time(it))];
    
    if (mod(it,15)==0||it==1)
        figure(1) % ----------------------------------------------------- %
        clf
        subplot(2,1,1)
        plotfield(D.T,M.X,M.Z,Pl,'contourf',...
            '\itT \rm\bf','quiver',ID.vx,ID.vz)
        subplot(2,1,2)
        plotfield(ID.v,M.X,M.Z,Pl,'pcolor',...
            '\itv \rm\bf')
    end
    
    %     % Capture the plot as an image
    %     frame       = getframe(h);
    %     im          = frame2im(frame);
    %     [imind,cm]  = rgb2ind(im,256);
    %
    %     % Write to the GIF File
    %     if it == 1
    %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    %     else
    %         imwrite(imind,cm,filename,'gif','WriteMode','append');
    %     end
    % ------------------------------------------------------------------- %
    
    %% Advektion ======================================================== %
    switch B.AdvMethod
        case 'semi-lag'
            if (it==1)
                %                 pause
                %                 keyboard
                [D.T]  =   SemiLagAdvection2D(ID,M,D.T,T.dt);
                
                % Speicher die alte Geschwindigkeit
                ID.vxo          =   ID.vx;
                ID.vzo          =   ID.vz;
            else
                [D.T]  =   SemiLagAdvection2D(ID,M,D.T,T.dt);
            end
    end
    % ------------------------------------------------------------------- %
    
    %% Diffusion ======================================================== %
    switch B.DiffMethod
        case 'explicit'
            D.T     =   SolveDiff2DexplicitSc(D.T,D.Q,T.dt,N,B);
        case 'implicit'
            D.T     =   SolveDiff2DimplicitSc(D.T,D.Q,T.dt,N,B);
            %         case 'implicito'
            %             [D.T,N] =   SolveDiff2Dimplicit_optSc(D.T,D.Q,T.dt,N,B);
        case 'ADI'
            D.T     =   SolveDiff2DADISc(D.T,D.Q,T.dt,N,B);
        case 'none'
        otherwise
            error('Diffusion scheme is not defined!')
    end
    % ------------------------------------------------------------------- %
    
    %% Wärmefluss an der Oberfläche ===================================== %
    D.dTtop         =   (D.T(1,:)-D.T(2,:))./N.dz;
    D.Nus(it)       =   trapz(M.x,D.dTtop);
    
    D.meanT(:,it)   =   mean(D.T,2);
    % ------------------------------------------------------------------- %
    
    %% Variable Viskositaet ============================================= %
    switch Py.eparam
        case 'variable'
            D.eta   =   exp( -Py.b.*((D.T-D.T(1,1))./(D.T(nz,1)-D.T(1,1)))...
                + Py.c.*M.Z);
    end
    % ------------------------------------------------------------------- %
    
    if (T.time(it) > T.tmax)
        disp('Maximale Zeit erreicht. Zeitschleife unterbrochen')
        T.indtime   = find(T.time(2:end)==0,1,'first');
        figure(1) % ----------------------------------------------------- %
        clf
        subplot(2,1,1)
        plotfield(D.T,M.X,M.Z,Pl,'contourf',...
            '\itT \rm\bf','quiver',ID.vx,ID.vz)
        subplot(2,1,2)
        plotfield(ID.v,M.X,M.Z,Pl,'pcolor',...
            '\itv \rm\bf')
        break
    elseif(it == T.itmax)
        T.indtime   = T.itmax;
        disp('Maximale Anzahl der Iterationen erreicht.')
        figure(1) % ----------------------------------------------------- %
        clf
        subplot(2,1,1)
        plotfield(D.T,M.X,M.Z,Pl,'contourf',...
            '\itT \rm\bf','quiver',ID.vx,ID.vz)
        subplot(2,1,2)
        plotfield(ID.v,M.X,M.Z,Pl,'pcolor',...
            '\itv \rm\bf')
        break
    end
end

%% Darstellen von Zeitparameters ======================================== %
figure(2)
subplot(3,1,1)
plot(T.time(1:T.indtime),D.Nus(1:T.indtime),...
    'LineWidth',2)
set(gca,'FontWeight','Bold');
xlabel('t'); ylabel('Nus')
title('Nusselt Number')

subplot(3,1,2)
plot(T.time(1:T.indtime),D.meanV(1:T.indtime),...
    'LineWidth',2)
set(gca,'FontWeight','Bold');
xlabel('t'); ylabel('VRMS')
title('Root Mean Square Velocity')

subplot(3,1,3)
plot(D.meanT(:,T.indtime),M.z,'LineWidth',2)
set(gca,'FontWeight','Bold');
xlabel('T'); ylabel('z')
title('Temperature Profile')
% saveas(figure(2),'TimeSeries','png')
% ----------------------------------------------------------------------- %
T.tend      = toc(T.tstart);

if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\DiffusionProblem')
    rmpath('..\AdvectionProblem')
    rmpath('..\StokesProblem')
    rmpath('..\SetUp')
else
    rmpath('../DiffusionProblem')
    rmpath('../AdvectionProblem')
    rmpath('../StokesProblem')
    rmpath('../SetUp')
end



% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


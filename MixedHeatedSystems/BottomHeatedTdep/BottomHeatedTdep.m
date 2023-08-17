clear
clc
% profile on

%% Add some paths ------------------------------------------------------- %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\DiffusionProblem')
    addpath('..\..\AdvectionProblem')
    addpath('..\..\StokesProblem')
    addpath('..\..\SetUp')
    addpath('..\..\ScaleParam')
else
    addpath('../../DiffusionProblem')
    addpath('../../AdvectionProblem')
    addpath('../../StokesProblem')
    addpath('../../SetUp')
    addpath('../../ScaleParam')
end

% T.tstart        =   tic;

%% Some initial definitions --------------------------------------------- %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
Py.Aparam       =   'temp';
Py.scale        =   'yes';

% Define method for solving the energy equation ------------------------- %
B.AdvMethod     =   'semi-lag';
B.DiffMethod    =   'explicit';
if ~exist('data/','dir')
    mkdir('data/')
end
% Variable oder konstante Vikositaet ------------------------------------ %
Py.eparam       =   'variable';
Py.b            =   log(1000); %log(16384);          % Temperaturabhaengigkeit
Py.c            =   0; %log(64);                  % Tiefenabhaengigkeit
B.EtaIni        =   'tdep';

% Define initial temperature anomaly ------------------------------------ %
B.Tini          =   'linano';
Py.tparam       =   'const';

% Define flow field ----------------------------------------------------- %
B.IniFlow       =   'none';
B.FlowFac       =   10;
% ----------------------------------------------------------------------- %

%% --------------------- Definition der Modelkonstanten ----------------- %
M.H         =   -400;          %   Modeltiefe [ in km ]
M.xmax      =   3;              %   Seitenverhaeltniss
% ----------------------------------------------------------------------- %

%% ------------------- Definition des Numerischen Gitters --------------- %
N.nz        =   101;             %   Vertikale Gitteraufloesung
N.nx        =   301;             %   Horizontale Gitteraufloesung
% ----------------------------------------------------------------------- %

%% -------------- Definition Physikalischer Konstanten ------------------ %
Py.g        =   9.81;               %   Schwerebeschleunigung [m/s^2]
Py.rho0     =   3300;               %   Hintergunddichte [kg/m^3]
Py.k        =   3;                  %   Thermische Leitfaehigkeit [ W/m/K ]
Py.cp       =   1000;               %   Heat capacity [ J/kg/K ]
Py.alpha    =   5e-5;               %   Thermischer Expnasionskoef. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]

Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]

Py.eta0     =   1e21;               %   Viskositaet [ Pa*s ]

Py.DeltaT   =   1000;           %   Temperaturdifferenz

% Falls Ra < 0 gesetzt ist, dann wird Ra aus den obigen Parametern
% berechnet. Falls Ra gegeben ist, dann wird die Referenzviskositaet so
% angepasst, dass die Skalierungsparameter die gegebene Rayleigh-Zahl
% ergeben.
Py.Ra       =   -999;           % Rayleigh number
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
T.tmaxini   =   4500;           %   Maximale Zeit in Ma
T.itmax     =   5000;           %   Maximal erlaubte Anzahl der Iterationen
T.dtfac     =   1.0;            %   Advektionscourantkriterium
T.dtdifac   =   0.9;            %   Diffusions Stabilitaetskriterium
% ----------------------------------------------------------------------- %
if Py.Ra < 0
    % Falls die Rayleigh Zahl nicht explizit angegeben wird, wird sie hier
    % berechnet.
    Py.Ra   =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H*1e3)^3/Py.eta0/Py.kappa;
else
    % Falls die Rayleigh Zahl gegeben ist m??ssen wir eine Variable
    % anpassen, z.B. die Referenzviskosit??t.
    Py.eta0 =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H*1e3)^3/Py.Ra/Py.kappa;
end
% ----------------------------------------------------------------------- %
%% ----------------------- Felddefinitionen ----------------------------- %
% Falls eine Strukturvariable schon zuvor definiert wurde, muss diese hier
% der Funktion als Eingabeargument mitgegeben werden. Das Selbe gilt fuer
% aller weiteren Funktionen.
[Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,N,M,T,Pl);
% ----------------------------------------------------------------------- %

%% ---------------- Definition der Anfangsbedingungen ------------------- %
[T,D,B,Ma,Py]           =   SetUpInitialConditions(T,D,Py,M,N,B);
% ----------------------------------------------------------------------- %
%% ------------------------- Plot Parameter ----------------------------- %
Pl.inc      =   min(N.nz/10,N.nx/5);
Pl.inc      =   round(Pl.inc);
Pl.xlab     =   '\bfx';
Pl.zlab     =   '\bfz';

switch Pl.plotfields
    case 'yes'
        if strcmp(getenv('OS'),'Windows_NT')
            set(figure(1),'position',[1.8,1.8,766.4,780.8]);
            h           =   figure(1);
        else
            set(figure(1),'position',[-1919,1,960,988]);
            h           =   figure(1);
        end
end

% Animation settings ---------------------------------------------------- %
switch Pl.savefig
    case 'yes'
        if ~exist('data/','dir')
            mkdir('data/')
        end
        Pl.filename    = ['data/ThermalConvect_Ra_',sprintf('%2.2e',Py.Ra),...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nz_',num2str(N.nz),'.gif'];
        % set(figure(1),'position',[55.4,125.8,1326.4,636.2]);
        h           =   figure(1);
end
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
    '\n Aufloesung (nx x nz): %i x %i',...
    '\n Refrenzviskositaet [Pa s]: %2.2e',...
    '\n  --------------------- ',...
    '\n '],Py.Ra,B.DiffMethod,B.AdvMethod,Py.eparam,B.Tini,B.IniFlow,...
    N.nx,N.nz,Py.eta0);
fprintf(['Maximum Time : %1.4g',...
    '\n Maximale Anzahl an Iterationen: %5i',...
    '\n  --------------------- \n'],...
    T.tmax,T.itmax)
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
    
    if (mod(it,10)==0||it==1)
        switch Pl.plotfields
            case 'yes'
                figure(1) % --------------------------------------------- %
                clf
                switch Py.eparam
                    case 'const'
                        ax1=subplot(2,1,1);
                        plotfield(D.T,M.X,M.Z,Pl,'contourf',...
                            '\itT \rm\bf','quiver',ID.vx,ID.vz)
                        colormap(ax1,flipud(Pl.lajolla))
                        caxis([0 1])
                        ax2=subplot(2,1,2);
                        plotfield(ID.v,M.X,M.Z,Pl,'pcolor',...
                            '\itv \rm\bf')
                        colormap(ax2,Pl.imola)
                    case 'variable'
                        ax1=subplot(2,1,1);
                        plotfield(D.T,M.X,M.Z,Pl,'contourf',...
                            '\itT\rm\bf','quiver',ID.vx,ID.vz)
                        colormap(ax1,flipud(Pl.lajolla))
                        ax2=subplot(2,1,2);
                        plotfield(log10(D.eta),M.X,M.Z,Pl,'contourf',...
                            '\it\eta\rm\bf','quiver',ID.vx,ID.vz)
                        colormap(ax2,flipud(Pl.lapaz))
                end
                
                switch Pl.savefig
                    case 'yes'
                        saveas(figure(1),...
                            ['data/ThermalConvect_Ra_',sprintf('%2.2e',Py.Ra),...
                            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
                            '_nz_',num2str(N.nz),'_',num2str(it)],'png')
                        
                        % Capture the plot as an image
                        frame       = getframe(h);
                        im          = frame2im(frame);
                        [imind,cm]  = rgb2ind(im,256);
                        
                        % Write to the GIF File
                        if it == 1
                            imwrite(imind,cm,Pl.filename,'gif', 'Loopcount',inf);
                        else
                            imwrite(imind,cm,Pl.filename,'gif','WriteMode','append');
                        end
                end
            case 'no'
                disp(['Iteration: ',sprintf('%i',it)])
        end
    end
    % ------------------------------------------------------------------- %
    
    %% Advektion ======================================================== %
    switch B.AdvMethod
        case 'semi-lag'
            if (it==1)
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
    
    %% W??rmefluss an der Oberfl??che ===================================== %
    D.dTtop         =   (D.T(1,:)-D.T(2,:))./N.dz;
    D.dTbot         =   (D.T(end-1,:)-D.T(end,:))./N.dz;
    %     D.Nus(it)       =   trapz(M.x,D.dTtop).*abs(M.H)/M.L;
    D.Nus(it)       =   mean(D.dTtop);
    
    D.meanT(:,it)   =   mean(D.T,2);
    % ------------------------------------------------------------------- %
    
    %% Variable Viskositaet ============================================= %
    switch B.EtaIni
        case 'tdep'
            D.eta   =   exp( -Py.b.*((D.T-D.T(1,1))./(D.T(N.nz,1)-D.T(1,1)))...
                + Py.c.*M.Z);
    end
    % ------------------------------------------------------------------- %
    [answer,T]  =   CheckBreakCriteria(it,T,D,M,Pl,ID,Py);
    switch answer
        case 'yes'
            break
    end
end
switch Pl.savefig
    case 'yes'
        saveas(figure(1),['data/ThermalConvect_Ra_',sprintf('%2.2e',Py.Ra),...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nz_',num2str(N.nz),'_SS'],'png')
end
%% Darstellen von Zeitparameters ======================================== %
PlotTimeSerieses(Py,T,D,M,N)
switch Pl.savefig
    case 'yes'
        saveas(figure(2),['data/ThermalConvect_Ra_',sprintf('%2.2e',Py.Ra),...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nz_',num2str(N.nz),'_TimeSeries'],'png')
end
% ----------------------------------------------------------------------- %
% T.tend      = toc(T.tstart);

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
switch Pl.savefig
    case 'yes'
        save('B.mat','B','-mat')
        save('D.mat','D','-mat')
        save('ID.mat','ID','-mat')
        save('M.mat','M','-mat')
        save('Ma.mat','Ma','-mat')
        save('N.mat','N','-mat')
        save('Pl.mat','Pl','-mat')
        save('Py.mat','Py','-mat')
        save('S.mat','S','-mat')
        save('T.mat','T','-mat')
end
% profile viewer

% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


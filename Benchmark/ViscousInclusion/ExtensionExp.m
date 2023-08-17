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

% eta2        = logspace(18,23,18);
% eta2        = fliplr(eta2); 
% 
% psiinc1     = zeros(length(eta2),1);
% psiinc2     = zeros(length(eta2),1);
% psiinc3     = zeros(length(eta2),1);
% 
% filename    = 'Inclusion_Deta_90.gif';
% 
% if strcmp(getenv('OS'),'Windows_NT')
%     set(figure(1),'position',[1.8,1.8,766.4,780.8]);
%     h           =   figure(1);
% else
%     set(figure(1),'position',[-1919,1,960,988]);
%     h           =   figure(1);
% end
% 
% for i = 1:length(eta2)
    
    %% Some initial definitions --------------------------------------------- %
    % Define method for solving the energy equation ------------------------- %
    B.AdvMethod     =   'none';
    B.DiffMethod    =   'none';
    
    % Variable oder konstante Vikositaet ------------------------------------ %
    Py.eparam       =   'variable';
    
    % Define initial temperature anomaly ------------------------------------ %
    B.Tini          =   'const';
    
    % Define flow field ----------------------------------------------------- %
    B.IniFlow       =   'xextension';
    B.FlowFac       =   [];
    % ----------------------------------------------------------------------- %
    
    %% --------------------- Definition der Modelkonstanten ----------------- %
    M.H         =   -1;           %   Modeltiefe [ in km ]
    M.xmax      =   1;               %   Seitenverhaeltniss
    % ----------------------------------------------------------------------- %
    
    %% ------------------- Definition des Numerischen Gitters --------------- %
    N.nz        =   51;             %   Vertikale Gitteraufloesung
    N.nx        =   51;             %   Horizontale Gitteraufloesung
    % ----------------------------------------------------------------------- %
    
    %% Tracer Advektionsmethode --------------------------------------------- %
    N.nmx       =   5;
    N.nmz       =   5;
    
    %% -------------- Definition Physikalischer Konstanten ------------------ %
    Py.g        =   10;                 %   Schwerebeschleunigung [m/s^2]
    Py.rho0     =   0;               %   Hintergunddichte [kg/m^3]
    Py.k        =   3;                  %   Thermische Leitfaehigkeit [ W/m/K ]
    Py.cp       =   1000;               %   Heat capacity [ J/kg/K ]
    Py.alpha    =   5e-5;               %   Thermischer Expnasionskoef. [ K^-1 ]
    
    Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]
    
    Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
    Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]
    
    Py.eta0     =   1e23;           % Viskositaet [ Pa*s ]
    Py.eta1     =   1e21;     % Inclusion viscosity
    
    Py.DeltaT   =   1000;           % Temperaturdifferenz
    % ----------------------------------------------------------------------- %
    
    %% -------------------- Definition der Randbedingungen ------------------ %
    % Geschwindigkeitsrandbedingungen --------------------------------------- %
    %   0 - no slip; 1 - free slip;
    B.tbc       =   1;              %   Obenw
    B.bbc       =   1;              %   Unten
    B.lbc       =   0;              %   Links
    B.rbc       =   0;              %   Rechts
    
    % Thermische Randbedingungen -------------------------------------------- %
    B.ttbc      =   'const';
    B.btbc      =   'const';
    B.ltbc      =   'const';
    B.rtbc      =   'const';
    
    % Waermerandbedinungen
    % Falls 'flux' - definiert es den Fluss
    % Falls 'const' (fuer top und bottom) - definiert es die Temperatur
    % Bottom, e.g. 1e-3
    B.lhf       =   1000;
    B.rhf       =   1000;
    B.thf       =   1000;
    B.bhf       =   1000;
    
    % Inklusionsbedingungen
    B.EtaIni        =   'ellipse';
    B.RotAng        =   90;             % positive -> counter clockwise
    B.EllA          =   -M.H/8*1e3;
    B.EllB          =   -M.H/4*1e3;
    % ----------------------------------------------------------------------- %
    
    %% ------------------------ Definition der Zeitkonstanten --------------- %
    T.tmax      =   4500;           %   Maximale Zeit in Ma
    T.itmax     =   1;            %   Maximal erlaubte Anzahl der Iterationen
    T.dtfac     =   1.0;            %   Advektionscourantkriterium
    T.dtdifac   =   1.0;            %   Diffusions Stabilitaetskriterium
    
    T.time      =   zeros(T.itmax,1);
    T.tmax      =   T.tmax.*1e6*(365.25*24*60*60);
    
    N.beenhere  =   0;
    N.beenhere2 =   0;
    % ----------------------------------------------------------------------- %
    
    %% ----------------------- Felddefinitionen ----------------------------- %
    % Falls eine Strukturvariable schon zuvor definiert wurde, muss diese hier
    % der Funktion als Eingabeargument mitgegeben werden. Das Selbe gilt fuer
    % aller weiteren Funktionen.
    [Py,D,ID,M,N,A]     =   SetUpFields(Py,N,M,T);
    % ----------------------------------------------------------------------- %
    
    %% ---------------- Definition der Anfangsbedingungen ------------------- %
    [D,B,Ma,Py] =   SetUpInitialConditions(D,Py,M,N,B);
    
    Py.Ra       =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.eta0/Py.kappa;
    % ----------------------------------------------------------------------- %
    
    %% ------------------------- Plot Parameter ----------------------------- %
    Pl.inc      =   min(N.nz/10,N.nx/10);
    Pl.inc      =   round(Pl.inc);
    % ----------------------------------------------------------------------- %
    
    % Information fuer das Befehlsfenster ----------------------------------- %
    %     fprintf([' Thermische Konvektion fuer Ra = %2.2e:\n  --------------------- ',...
    %         '\n Diffusion mit: %s',...
    %         '\n Advektion mit: %s',...
    %         '\n Viskositaet ist: %s',...
    %         '\n Anfangstemperaturfeld: %s',...
    %         '\n Anfangsgeschwindigkeitsfeld: %s',...
    %         '\n  --------------------- ',...
    %         '\n '],Py.Ra,B.DiffMethod,B.AdvMethod,Py.eparam,B.Tini,B.IniFlow);
    %     fprintf(['Maximum Time : %1.4g [ Ma ]',...
    %         '\n Maximale Anzahl an Iterationen: %5i',...
    %         '\n  --------------------- \n'],...
    %         T.tmax/(1e6*(365.25*24*60*60)),T.itmax)
    %     % pause
    
    % ----------------------------------------------------------------------- %
    
    %% BEGINN DER ZEITSCHLEIFE ============================================= %%
    for it = 1:T.itmax
        %% Erstellung der Koeffizienten Matrix und rechten Seite des
        %     if(strcmp(B.AdvMethod,'none')==0)
        switch Py.eparam
            case 'const'
                [D,A]       =   solveSECE_const_Eta(D,Py,N,B,A);
                if (it == 2)
                    N.beenhere = 1;
                end
            case 'variable'
                [D,A]       =   solveSECE(D,Py,N,B);
            otherwise
                error('Viskositaet nicht definiert! Siehe Parameter Py.eparam.')
        end
        %     end
        % ------------------------------------------------------------------- %
        
        %% Berechnung der Zeitschrittlaenge --------------------------------- %
        T.dt        =   T.dtfac*min(N.dx,abs(N.dz))/...
            (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));
        
        T.dtdiff    =   T.dtdifac*min([N.dx,abs(N.dz)])^2/Py.kappa/4;
        
        T.dt        =   min(T.dt,T.dtdiff);
        
        if it>1
            T.time(it)  =   T.time(it-1) + T.dt;
        end
        % ------------------------------------------------------------------- %
        
        %% Interpolation der Geschwindigkeiten auf das regulaere Gitter ----- %
        [ID]        =   InterpStaggered(D,ID,N,'velocity');
        D.meanV(it) = mean(ID.v,'all');   % Mittleregeschwindigkeit
        % ------------------------------------------------------------------- %
        
        [ID]        =   GetStrainRate(ID,N);
        
        ID.tauII    =   ID.eII.*Py.eta.*2;
        
        ID.psi      =   ID.eII.*ID.tauII;
        
        incind      =   log10(Py.eta)==log10(Py.eta1);                
        
        %% Darstellung der Daten -------------------------------------------- %
        Pl.time     =   '';
        
        if (mod(it,5)==0||it==1)
            figure(1) % ----------------------------------------------------- %
            clf
            subplot(2,2,1)
            plotfield(log10(Py.eta),M.X./1000,M.Z./1000,Pl,'pcolor',...
                '\itlog_{10} ( \eta ) [Pa s]\rm\bf [ K ]','quiver',ID.vx,ID.vz)
            subplot(2,2,2)
            plotfield(log10(ID.psi),M.X./1000,M.Z./1000,Pl,'pcolor',...
                '\it log_{10} ( \psi ) [Pa/s]\rm\bf')
            subplot(2,2,3)
            plotfield(log10(ID.eII),M.X./1000,M.Z./1000,Pl,'pcolor',...
                '\itlog_{10} ( \epsilon_{II} ) [s^{-1}]\rm\bf')
            subplot(2,2,4)
            plotfield(ID.tauII./1e6,M.X./1000,M.Z./1000,Pl,'pcolor',...
                '\it \tau_{II} [MPa]\rm\bf')
        end
        
%         % Capture the plot as an image
%         frame       = getframe(h);
%         im          = frame2im(frame);
%         [imind,cm]  = rgb2ind(im,256);
%         
%         % Write to the GIF File
%         if i == 1
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append');
%         end
        % ------------------------------------------------------------------- %
        
    end
    
%     psiinc1(i) = mean(ID.psi(incind));
%     psiinc2(i) = -sum(ID.psi(incind))./sum(incind(:).*N.dx.*N.dz);
%     psiinc3(i) = sum(ID.psi(incind))./(pi/B.EllA/B.EllB);
%     
% end
% 
% set(figure(2),'position',[488,162.6,669.8,599.4000])
% figure(2)
% clf
% plot(eta2./Py.eta0,psiinc1,'LineWidth',2)
% hold on
% plot(eta2./Py.eta0,psiinc2,'LineWidth',2)
% plot(eta2./Py.eta0,psiinc3,'LineWidth',2)
% xlabel('\eta_{inc}/\eta_{out}'); ylabel('\psi')
% legend('\langle\psi\rangle','\psi_{integral1}','\psi_{integral2}',...
%     'Location','SouthEast')
% set(gca,'FontWeight','Bold','xscale','log','yscale','log','LineWidth',2,'FontSize',15)
% % axis([0 180 1e-11 1e-6])
% 
% % saveas(figure(2),'Dissipation_Deta_90','png')

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


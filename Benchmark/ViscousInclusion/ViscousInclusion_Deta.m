clear
clc

%% Add some paths ------------------------------------------------------- %
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

angle = [0 22.5 45 90];

for k = 1:length(angle)
    
    
    eta2        = logspace(18,28,40);
    eta2        = fliplr(eta2);
    
    psiinc1ma   = zeros(length(eta2),1);
    eIIincma    = zeros(length(eta2),1);
    tauIIincma  = zeros(length(eta2),1);
    
    psiinc1mi   = zeros(length(eta2),1);
    eIIincmi    = zeros(length(eta2),1);
    tauIIincmi  = zeros(length(eta2),1);
    
    psiinc1std  = zeros(length(eta2),1);
    eIIincstd   = zeros(length(eta2),1);
    tauIIincstd = zeros(length(eta2),1);
    
    psiinc1     = zeros(length(eta2),1);
    eIIinc      = zeros(length(eta2),1);
    tauIIinc    = zeros(length(eta2),1);
    
    psimat1ma   = zeros(length(eta2),1);
    eIImatma    = zeros(length(eta2),1);
    tauIImatma  = zeros(length(eta2),1);
    
    psimat1mi   = zeros(length(eta2),1);
    eIImatmi    = zeros(length(eta2),1);
    tauIImatmi  = zeros(length(eta2),1);
    
    psimat1std  = zeros(length(eta2),1);
    eIImatstd   = zeros(length(eta2),1);
    tauIImatstd = zeros(length(eta2),1);
    
    psimat1     = zeros(length(eta2),1);
    eIImat      = zeros(length(eta2),1);
    tauIImat    = zeros(length(eta2),1);
    
    Orientation = angle(k);
    
    if strcmp(getenv('OS'),'Windows_NT')
        set(figure(1),'position',[1.8,1.8,766.4,780.8]);
        h           =   figure(1);
    else
        set(figure(1),'position',[-1919,1,960,988]);
        h           =   figure(1);
    end
    
    Pl.savefig      =   'no';
    Pl.plotfields   =   'yes';
    
    for i = 1:length(eta2)
        
        %% Some initial definitions ------------------------------------- %
        % Define method for solving the energy equation ----------------- %
        B.AdvMethod     =   'none';
        B.DiffMethod    =   'none';
        B.Aparam        =   'none';
        Py.scale        =   'none';
        
        % Variable oder konstante Vikositaet ---------------------------- %
        Py.eparam       =   'variable';
        Py.tparam       =   '';
        
        % Define initial temperature anomaly ---------------------------- %
        B.Tini          =   'const';
        
        % Define flow field --------------------------------------------- %
        B.IniFlow       =   'PureShear';
        B.FlowFac       =   [];
        % --------------------------------------------------------------- %
        
        %% --------------------- Definition der Modelkonstanten --------- %
        M.H         =   -1;             %   Modeltiefe [ in km ]
        M.xmax      =   1;              %   Seitenverhaeltniss
        % --------------------------------------------------------------- %
        
        %% ------------------- Definition des Numerischen Gitters ------- %
        N.nz        =   51;            %   Vertikale Gitteraufloesung
        N.nx        =   51;            %   Horizontale Gitteraufloesung
        % --------------------------------------------------------------- %
        
        %% Tracer Advektionsmethode ------------------------------------- %
        N.nmx       =   5;
        N.nmz       =   5;
        
        %% -------------- Definition Physikalischer Konstanten ---------- %
        Py.g        =   10;             %   Schwerebeschleunigung [m/s^2]
        Py.rho0     =   3200;           %   Hintergunddichte [kg/m^3]
        Py.k        =   3;              %   Thermische Leitfaehigkeit [ W/m/K ]
        Py.cp       =   1000;           %   Heat capacity [ J/kg/K ]
        Py.alpha    =   5e-5;           %   Thermischer Expnasionskoef. [ K^-1 ]
        
        Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]
        
        Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
        Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]
        
        Py.eta0     =   1e23;           % Viskositaet [ Pa*s ]
        Py.eta1     =   eta2(i);        % Inclusion viscosity
        Py.rho1     =   Py.rho0; 
        
        Py.DeltaT   =   1000;           % Temperaturdifferenz
        % --------------------------------------------------------------- %
        
        %% -------------------- Definition der Randbedingungen ---------- %
        % Geschwindigkeitsrandbedingungen ------------------------------- %
        %   0 - no slip; 1 - free slip;
        B.tbc       =   0;              %   Obenw
        B.bbc       =   0;              %   Unten
        B.lbc       =   0;              %   Links
        B.rbc       =   0;              %   Rechts
        
        % Thermische Randbedingungen ------------------------------------ %
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
        B.ebg           =   -1e-15;         % < 0 compression
        B.RotAng        =   Orientation;    % positive -> counter clockwise
        B.EllA          =   3e2;            % [ m ]
        B.EllB          =   1e2;          % [ m ]
        B.T0            =   1000; 
        B.TAmpl         =   1000; 
        
        switch Pl.savefig
            case 'yes'
                ModDir      = ['data/Ell_a_',num2str(B.EllA),'_b_',...
                    num2str(B.EllB)];
                if ~exist(ModDir,'dir')
                    mkdir(ModDir)
                end
                
                filename    = [ModDir,'/Inclusion_Deta_',B.IniFlow,'_',...
                    num2str(Orientation),'.gif'];
        end
        % --------------------------------------------------------------- %
        
        %% ------------------------ Definition der Zeitkonstanten ------- %
        T.tmaxini   =   4500;           %   Maximale Zeit in Ma
        T.itmax     =   1;              %   Maximal erlaubte Anzahl der Iterationen
        T.dtfac     =   1.0;            %   Advektionscourantkriterium
        T.dtdifac   =   1.0;            %   Diffusions Stabilitaetskriterium
        % --------------------------------------------------------------- %
        
        %% ----------------------- Felddefinitionen --------------------- %
        % Falls eine Strukturvariable schon zuvor definiert wurde, muss diese hier
        % der Funktion als Eingabeargument mitgegeben werden. Das Selbe gilt fuer
        % aller weiteren Funktionen.
        [Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
        % --------------------------------------------------------------- %
        
        %% ---------------- Definition der Anfangsbedingungen ----------- %
        [T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
        
        Py.Ra       =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.eta0/Py.kappa;
        % --------------------------------------------------------------- %
        
        %% ------------------------- Plot Parameter --------------------- %
        Pl.inc      =   min(N.nz/10,N.nx/10);
        Pl.inc      =   round(Pl.inc);
        % --------------------------------------------------------------- %        
        %% Scale Parameters ============================================= %
        [M,N,D,T,S]         =   ScaleParameters(B,M,Py,N,D,T);
        switch B.IniFlow
            case 'SimpleShear'
                gr  =   1;
                er  =   0;
            case 'PureShear'
                gr  =   0;
                er  =   -1;
        end
        
        [ Vx_N,Vx_S,Vx_W,Vx_E,Vz_N,Vz_S,Vz_W,Vz_E, Pa, Vxa, Vza ] = ...
            Dani_Solution_vec(M.x-M.L/2,M.z-M.H/2,M.x1-M.L/2,M.z1-M.H/2,...
            (B.EllA)/(-M.H*1e3),Py.eta1/Py.eta0,N.nx1,N.nz1,...
            gr,er);
        D.Pa =    Pa'; D.Vxa =  Vxa'; D.Vza =   Vza';
        clear Vxa Vza Pa
        % Overwrite BC
        Vx_S(1)     = Vx_W(1);  Vx_S(end)     = Vx_E(1);
        Vx_N(1)     = Vx_W(end);Vx_N(end)     = Vx_E(end);
        D.vxi(1:end-1,1)    =   Vx_W; D.vxi(1:end-1,N.nx) =   Vx_E;
        D.vxi(1,:)          =   Vx_S; D.vxi(N.nz1,:)      =   Vx_N;
        D.vzi(1:end,1)      =   Vz_W; D.vzi(1:end,N.nx1)  =   Vz_E;
        D.vzi(1,1:end-1)    =   Vz_S; D.vzi(N.nz,1:end-1) =   Vz_N;
        % keyboard
        
        
        %% BEGINN DER ZEITSCHLEIFE ===================================== %%
        for it = 1:T.itmax
            %% Erstellung der Koeffizienten Matrix und rechten Seite des
            %     if(strcmp(B.AdvMethod,'none')==0)
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
            %     end
            % ----------------------------------------------------------- %
            
            %% Interpolation der Geschwindigkeiten auf das regulaere Gitter - %
            [ID]        =   InterpStaggered(D,ID,N,'velocity');
            D.meanV(it) = mean(ID.v,'all');   % Mittleregeschwindigkeit
            % ----------------------------------------------------------- %
            [ID]        =   GetStrainRate(ID,N);
            ID.tauII    =   ID.eII.*D.eta.*2;
            ID.psi      =   ID.eII.*ID.tauII;
            incind      =   log10(D.eta)==log10(Py.eta1/Py.eta0);
            matind      =   log10(D.eta)==log10(Py.eta0/Py.eta0);
            
            %% Darstellung der Daten ------------------------------------ %
            Pl.time     =   '';
            Pl.xlab     =   '$$x$$';
            Pl.zlab     =   '$$z$$';
            
            if (mod(it,5)==0||it==1)
                figure(1) % --------------------------------------------- %
                clf
                ax1 = subplot(2,2,1);
                plotfield(log10(D.eta),M.X,M.Z,Pl,'pcolor',...
                    '$$log_{10}\ (\eta)$$','quiver',ID.vx,ID.vz)
                colormap(ax1,Pl.lapaz)
                ax2 = subplot(2,2,2);
                plotfield(log10(ID.psi),M.X,M.Z,Pl,'pcolor',...
                    '$$log_{10}\ (\psi)$$')
                colormap(ax2,Pl.imola)
                ax3 = subplot(2,2,3);
                plotfield(log10(ID.eII),M.X,M.Z,Pl,'pcolor',...
                    '$$log_{10}\ ( \varepsilon_{II} )$$')
                colormap(ax3,Pl.batlowW)
                ax4 = subplot(2,2,4);
                plotfield(log10(ID.tauII),M.X,M.Z,Pl,'pcolor',...
                    '$$log_{10}\ ( \tau_{II} )$$')
                colormap(ax4,Pl.nuuk)
            end
            
            switch Pl.savefig
                case 'yes'
                    % Capture the plot as an image
                    frame       = getframe(h);
                    im          = frame2im(frame);
                    [imind,cm]  = rgb2ind(im,256);
                    
                    % Write to the GIF File
                    if i == 1
                        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                    else
                        imwrite(imind,cm,filename,'gif','WriteMode','append');
                    end
            end
            % ----------------------------------------------------------- %
            
        end
        
        psiinc1(i)      = mean(ID.psi(incind));
        eIIinc(i)       = mean(ID.eII(incind));
        tauIIinc(i)     = mean(ID.tauII(incind));
        
        psiinc1ma(i)    = max(ID.psi(incind));
        eIIincma(i)     = max(ID.eII(incind));
        tauIIincma(i)   = max(ID.tauII(incind));
        
        psiinc1mi(i)    = min(ID.psi(incind));
        eIIincmi(i)     = min(ID.eII(incind));
        tauIIincmi(i)   = min(ID.tauII(incind));
        
        psiinc1std(i)   = std(ID.psi(incind));
        eIIincstd(i)    = std(ID.eII(incind));
        tauIIincstd(i)  = std(ID.tauII(incind));
        
        % -------
        psimat1(i)      = mean(ID.psi(matind));
        eIImat(i)       = mean(ID.eII(matind));
        tauIImat(i)     = mean(ID.tauII(matind));
        
        psimat1ma(i)    = max(ID.psi(matind));
        eIImatma(i)     = max(ID.eII(matind));
        tauIImatma(i)   = max(ID.tauII(matind));
        
        psimat1mi(i)    = min(ID.psi(matind));
        eIImatmi(i)     = min(ID.eII(matind));
        tauIImatmi(i)   = min(ID.tauII(matind));
        
        psimat1std(i)   = std(ID.psi(matind));
        eIImatstd(i)    = std(ID.eII(matind));
        tauIImatstd(i)  = std(ID.tauII(matind));
        
    end
    data    =   [eta2'./Py.eta0 ...
        psiinc1 psiinc1ma psiinc1mi psiinc1std ...
        eIIinc eIIincma eIIincmi eIIincstd ...
        tauIIinc tauIIincma tauIIincmi tauIIincstd ...
        psimat1 psimat1ma psimat1mi psimat1std ...
        eIImat eIImatma eIImatmi eIImatstd ...
        tauIImat tauIImatma tauIImatmi tauIImatstd  ...
        ];
    switch Pl.savefig
        case 'yes'
            save([ModDir,'/Data_',B.IniFlow,'_',...
                num2str(Orientation),'.mat'],'data','-mat')
    end
    
    set(figure(2),'position',[488,162.6,669.8,599.4000])
    figure(2)
    clf
    yyaxis left; plot(eta2./Py.eta0,psiinc1,'LineWidth',2)
    hold on
    plot(eta2./Py.eta0,tauIIinc,'LineWidth',2)
    ylabel('$$\psi or \tau_{II}$$','Interpreter','latex')
    set(gca,'FontWeight','Bold','xscale','log','yscale','log',...
        'LineWidth',2,'FontSize',15,'TickLabelInterpreter','latex')
    yyaxis right; plot(eta2./Py.eta0,eIIinc,'LineWidth',2)
    set(gca,'FontWeight','Bold','xscale','log','yscale','log',...
        'LineWidth',2,'FontSize',15,'TickLabelInterpreter','latex')
    ylabel('$$\varepsilon_{II}$$','Interpreter','latex')
    xlabel('$$\eta_{inc}/\eta_{out}$$','Interpreter','latex');
    legend('$$\langle \psi \rangle$$','$$\langle \tau_{II} \rangle$$',...
        '$$\langle \epsilon_{II} \rangle$$','Location','SouthEast',...
        'Interpreter','latex')
    % axis([0 180 1e-11 1e-6])
    switch Pl.savefig
        case 'yes'
            saveas(figure(2),[ModDir,'/Dissipation_Deta_',B.IniFlow,'_',...
                num2str(Orientation)],'png')
    end
    
    set(figure(3),'Position',[1.8000,1.8000,766.4000,780.8000])
    
    figure(3)
    subplot(3,1,3)
    plot(eta2./Py.eta0,psiinc1,'LineWidth',2)
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \psi_{II} \rangle$$','Interpreter','latex')
    set(gca,'LineWidth',2,'FontSize',20,...
        'yscale','log','xscale','log','TickLabelInterpreter','latex');
    % axis([1e-4 1e2 1e-16 1e-11])
    subplot(3,1,2)
    plot(eta2./Py.eta0,tauIIinc,'LineWidth',2)
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \tau_{II} \rangle$$','Interpreter','latex')
    set(gca,'LineWidth',2,'FontSize',20,...
        'yscale','log','xscale','log','TickLabelInterpreter','latex');
    subplot(3,1,1)
    plot(eta2./Py.eta0,eIIinc,'LineWidth',2)
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \epsilon_{II} \rangle $$','Interpreter','latex')
    set(gca,'LineWidth',2,'FontSize',20,...
        'yscale','log','xscale','log','TickLabelInterpreter','latex');
    switch Pl.savefig
        case 'yes'
            saveas(figure(3),[ModDir,'/eps_tau_psi_Deta_',B.IniFlow,'_',...
                num2str(Orientation)],'png')
    end
    
    figure(4)
    plot(eIIinc,tauIIinc,'LineWidth',2)
    xlabel('$$\langle \epsilon_{II} \rangle$$','Interpreter','latex')
    ylabel('$$\langle \tau_{II} \rangle$$','Interpreter','latex')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,...
        'yscale','log','xscale','log','TickLabelInterpreter','latex');
    % axis([1e-16 1e-10 1e-2 1e3]);
    axis square
    switch Pl.savefig
        case 'yes'
            saveas(figure(4),[ModDir,'/eps_tau',B.IniFlow,'_',...
                num2str(Orientation)],'png')
    end
    close all
end

if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\..\AdvectionProblem')
    rmpath('..\..\StokesProblem')
    rmpath('..\..\SetUp')
    rmpath('..\..\ScaleParam')
else
    rmpath('../../AdvectionProblem')
    rmpath('../../StokesProblem')
    rmpath('../../SetUp')
    rmpath('../../ScaleParam')
end

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

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

rotAngle    = 0:2:180;
eta2        = 1e28;

psiinc1     = zeros(length(rotAngle),1);
eIIinc      = zeros(length(rotAngle),1);
tauIIinc    = zeros(length(rotAngle),1);

if strcmp(getenv('OS'),'Windows_NT')
    set(figure(1),'position',[1.8,1.8,766.4,780.8]);
    h           =   figure(1);
else
    set(figure(1),'position',[-1919,1,960,988]);
    h           =   figure(1);
end

Pl.savefig      =   'no';
Pl.plotfields   =   'yes';

for i = 1:length(rotAngle)
    
    %% Some initial definitions ----------------------------------------- %
    % Define method for solving the energy equation --------------------- %
    B.AdvMethod     =   'none';
    B.DiffMethod    =   'none';
    B.Aparam       =   'none';
    Py.scale        =   'yes';
    
    % Variable oder konstante Vikositaet -------------------------------- %
    Py.eparam       =   'variable';
    Py.tparam       =   '';
    
    % Define initial temperature anomaly -------------------------------- %
    B.Tini          =   'const';
    
    % Define flow field ------------------------------------------------- %
    B.IniFlow       =   'SimpleShear';
    B.FlowFac       =   [];
    % ------------------------------------------------------------------- %
    
    %% --------------------- Definition der Modelkonstanten ------------- %
    M.H         =   -1;             %   Modeltiefe [ in km ]
    M.xmax      =   1;              %   Seitenverhaeltniss
    % ------------------------------------------------------------------- %
    
    %% ------------------- Definition des Numerischen Gitters ----------- %
    N.nz        =   201;            %   Vertikale Gitteraufloesung
    N.nx        =   201;            %   Horizontale Gitteraufloesung
    % ------------------------------------------------------------------- %
    
    %% Tracer Advektionsmethode ----------------------------------------- %
    N.nmx       =   5;
    N.nmz       =   5;
    
    %% -------------- Definition Physikalischer Konstanten -------------- %
    Py.g        =   10;                 %   Schwerebeschleunigung [m/s^2]
    Py.rho0     =   3200;               %   Hintergunddichte [kg/m^3]
    Py.k        =   3;                  %   Thermische Leitfaehigkeit [ W/m/K ]
    Py.cp       =   1000;               %   Heat capacity [ J/kg/K ]
    Py.alpha    =   5e-5;               %   Thermischer Expnasionskoef. [ K^-1 ]
    
    Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]
    
    Py.Q0       =   0;              % Waermeproduktionsrate pro Volumen [W/m^3]
    Py.Q0       =   Py.Q0/Py.rho0;  % Waermeproduktionsrate pro Masse [W/kg]
    
    Py.eta0     =   1e23;           % Viskositaet [ Pa*s ]
    Py.eta1     =   eta2;           % Inclusion viscosity
    
    Py.DeltaT   =   1000;           % Temperaturdifferenz
    % ------------------------------------------------------------------- %
    
    %% -------------------- Definition der Randbedingungen -------------- %
    % Geschwindigkeitsrandbedingungen ----------------------------------- %
    %   0 - no slip; 1 - free slip;
    B.tbc       =   0;              %   Obenw
    B.bbc       =   0;              %   Unten
    B.lbc       =   0;              %   Links
    B.rbc       =   0;              %   Rechts
    
    % Thermische Randbedingungen ---------------------------------------- %
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
    B.ebg           =   -1e-15;         % < 0 compression; is overwriten by analytical solution!
    B.RotAng        =   rotAngle(i);             % positive -> counter clockwise
    B.EllA          =   2e2;            % [ m ]
    B.EllB          =   1e2;
    % ------------------------------------------------------------------- %
    
    %% ------------------------ Definition der Zeitkonstanten ----------- %
    T.tmaxini   =   4500;           %   Maximale Zeit in Ma
    T.itmax     =   1;              %   Maximal erlaubte Anzahl der Iterationen
    T.dtfac     =   1.0;            %   Advektionscourantkriterium
    T.dtdifac   =   1.0;            %   Diffusions Stabilitaetskriterium
    % ------------------------------------------------------------------- %
    
    %% ----------------------- Felddefinitionen ------------------------- %
    % Falls eine Strukturvariable schon zuvor definiert wurde, muss diese hier
    % der Funktion als Eingabeargument mitgegeben werden. Das Selbe gilt fuer
    % aller weiteren Funktionen.
    [Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
    % ------------------------------------------------------------------- %
    
    %% ---------------- Definition der Anfangsbedingungen --------------- %
    [T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
    
    Py.Ra       =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.eta0/Py.kappa;
    % ------------------------------------------------------------------- %
    
    %% ------------------------- Plot Parameter ------------------------- %
    Pl.inc      =   min(N.nz/10,N.nx/10);
    Pl.inc      =   round(Pl.inc);
    filename    = ['data/Inclusion_Rotation_m_',num2str(Py.eta1/Py.eta0),...
        B.IniFlow,'.gif'];
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    
    %% Scale Parameters ================================================= %
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
    
    
    %% BEGINN DER ZEITSCHLEIFE ========================================= %%
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
        % --------------------------------------------------------------- %
        
        %% Interpolation der Geschwindigkeiten auf das regulaere Gitter - %
        [ID]        =   InterpStaggered(D,ID,N,'velocity');
        D.meanV(it) = mean(ID.v,'all');   % Mittleregeschwindigkeit
        % --------------------------------------------------------------- %
        [ID]        =   GetStrainRate(ID,N);
        ID.tauII    =   ID.eII.*D.eta.*2;
        ID.psi      =   ID.eII.*ID.tauII;
        incind      =   log10(D.eta)==log10(Py.eta1/Py.eta0);
        
        %% Darstellung der Daten ---------------------------------------- %
        Pl.time     =   '';
        Pl.xlab     =   'x';
        Pl.zlab     =   'z';
        
        if (mod(it,5)==0||it==1)
            figure(1) % ------------------------------------------------- %
            clf
            subplot(2,2,1)
            plotfield(log10(D.eta),M.X,M.Z,Pl,'pcolor',...
                '\itlog_{10} ( \eta ) \rm\bf','quiver',ID.vx,ID.vz)
            subplot(2,2,2)
            plotfield(log10(ID.psi),M.X,M.Z,Pl,'pcolor',...
                '\it log_{10} ( \psi ) \rm\bf')
            subplot(2,2,3)
            plotfield(log10(ID.eII),M.X,M.Z,Pl,'pcolor',...
                '\itlog_{10} ( \epsilon_{II} ) \rm\bf')
            subplot(2,2,4)
            plotfield(log10(ID.tauII),M.X,M.Z,Pl,'pcolor',...
                '\it log_{10} ( \tau_{II} ) \rm\bf')
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
        % --------------------------------------------------------------- %
        
    end
    psiinc1(i) = mean(ID.psi(incind));
    eIIinc(i)   = mean(ID.eII(incind));
    tauIIinc(i) = mean(ID.tauII(incind));
    
end

set(figure(2),'position',[488,162.6,669.8,599.4000])
figure(2)
clf
yyaxis left; plot(rotAngle,psiinc1,'LineWidth',2)
hold on
plot(rotAngle,tauIIinc,'LineWidth',2)
xlabel('\alpha'); ylabel('\psi or \tau_{II}')
set(gca,'FontWeight','Bold','yscale','log','LineWidth',2,'FontSize',15)
yyaxis right; plot(rotAngle,eIIinc,'LineWidth',2)
set(gca,'FontWeight','Bold','yscale','log','LineWidth',2,'FontSize',15)
ylabel('\epsilon_{II}'); xlabel('\alpha');
legend('\langle\psi\rangle','\langle\tau_{II}\rangle',...
    '\langle\epsilon_{II}\rangle','Location','SouthEast')
% axis([0 180 1e-6 5e-5])

switch Pl.savefig
    case 'yes'
        saveas(figure(2),['data/Dissipation_Rotation',B.IniFlow],'png')
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

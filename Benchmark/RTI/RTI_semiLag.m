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
else
    addpath('../../DiffusionProblem')
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
Py.scale        =   'no';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'semi-lag';
B.Aparam        =   'comp';
B.DiffMethod    =   'none';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'variable';
B.EtaIni        =   'RTI';
B.lambda 	=   0.5; 		% Wellenlaenge der Perturbation [ km ]
B.deltaA        =   100;        % Amplitude [ km ]
% ======================================================================= %
%% ================== Define initial temperature anomaly ================ %
B.Tini          =   'none';
Py.tparam       =   'none';
% ======================================================================= %
%% ========================= Define flow field ========================== %
B.IniFlow       =   'none';
B.FlowFac       =   10;
% ======================================================================= %
%% ==================== Define model geometry constants ================= %
M.H         =   -1;                 %   Modeltiefe [ in km ]
M.xmax      =   -2*B.lambda/M.H;   %   Seitenverhaeltniss
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   101;                 %   Vertikale Gitteraufloesung
N.nx        =   (N.nz-1)*M.xmax+1;   %   Horizontale Gitteraufloesung
% ======================================================================= %
%% ====================== Tracer advection method ======================= %
% So far only required for initialization!
N.nmx       =   5;
N.nmz       =   5;
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   10.0;           %   Schwerebeschleunigung [m/s^2]
Py.rho0     =   3300;           %   Hintergunddichte [kg/m^3]
Py.k        =   3;              %   Thermische Leitfaehigkeit [ W/m/K ]
Py.cp       =   1000;           %   Heat capacity [ J/kg/K ]
Py.alpha    =   5e-5;           %   Thermischer Expnasionskoef. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]

Py.Q0       =   0;              %   Waermeproduktionsrate pro Volumen [W/m^3]

Py.eta0     =   1e21;           %   Viskositaet comp. 0 [ Pa*s ]
Py.eta1     =   1e18;           %   Viskositaet comp. 1 [ Pa*s ]

Py.rho1     =   3000;           %   Density comp. 1 [ kg/m^3 ]
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Geschwindigkeitsrandbedingungen --------------------------------------- %
%   0 - no slip; 1 - free slip;
B.tbc       =   0;          %   Top
B.bbc       =   0;          %   Bottom
B.lbc       =   1;          %   Left
B.rbc       =   1;          %   Right
% ======================================================================= %
%% ====================== Define time constants ========================= %
T.tmaxini   =   4500;       %   Maximum time in Ma
T.itmax     =   50;         %   Maximal erlaubte Anzahl der Iterationen
T.dtfac     =   1.0;        %   Advektionscourantkriterium
% ======================================================================= %
%% ========================= Define fields required ===================== %
[Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
% ======================================================================= %
%% ======================== Setup initial conditions ==================== %
[T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
% ======================================================================= %
%% ========================= Plot parameter ============================= %
Pl.inc      =   min(N.nz/10,N.nx/5);
Pl.inc      =   round(Pl.inc);
Pl.xlab     =   '\bfx [ km ]';
Pl.zlab     =   '\bfz [ km ]';
switch Pl.plotfields
    case 'yes'
        if strcmp(getenv('OS'),'Windows_NT')
            set(figure(3),'position',[1.8,1.8,766.4,780.8]);
            h           =   figure(3);
        else
            set(figure(3),'position',[-1919,1,960,988]);
            h           =   figure(3);
        end
end
% Animation settings ---------------------------------------------------- %
switch Pl.savefig
    case 'yes'
        M.ModDir    =   ['data/RTI_adv_',B.AdvMethod,...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nz_',num2str(N.nz),'/'];
        if ~exist(M.ModDir,'dir')
            mkdir(M.ModDir)
        end
        Pl.filename    = [M.ModDir,'Evolution.gif'];
        h           =   figure(3);
end
% ======================================================================= %
%% ========================== Scale Parameters ========================== %
switch lower(Py.scale)
    case 'yes'
        [M,N,D,T,S]         =   ScaleParameters(B,M,Py,N,D,T);
end
% ======================================================================= %
%% ================ Information for the command window ================== %
fprintf([' Rayleigh-Taylor Instabilitaet  --------------------- ',...
    '\n Advektion mit: %s',...
    '\n Viskositaet ist: %s',...
    '\n Anfangstemperaturfeld: %s',...
    '\n Anfangsgeschwindigkeitsfeld: %s',...
    '\n Aufloesung (nx x nz): %i x %i',...
    '\n Refrenzviskositaet [Pa s]: %2.2e',...
    '\n  --------------------- ',...
    '\n\n '],B.AdvMethod,Py.eparam,B.Tini,B.IniFlow,...
    N.nx,N.nz,Py.eta0);
fprintf(['Maximum Time : %1.4g',...
    '\n Maximale Anzahl an Iterationen: %5i',...
    '\n  --------------------- \n'],...
    T.tmax/1e6/(365.25*24*60*60),T.itmax)
% ======================================================================= %
%% ========================= Time loop ================================= %%
for it = 1:T.itmax
    %     disp(['Iteration: ',sprintf('%i',it)])
    if(strcmp(B.AdvMethod,'none')==0)
        switch Py.eparam
            case 'const'
                [D,A]       =   solveSECE_const_Eta(D,Py,N,B,A);
                if (it == 2)
                    N.beenhere = 1;
                end
            case 'variable'
                [D,A]       =   solveSECE(D,Py,N,B);
            otherwise
                error('Viscosity not difined! Check Py.eparam parameter.')
        end
    end
    % =================================================================== %
    %% ===================== Calculate time stepping ==================== %
    T.dt        =   T.dtfac*min(N.dx,abs(N.dz))/...
        (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));    
    if it>1
        T.time(it)  =   T.time(it-1) + T.dt;
    end
    % =================================================================== %
    %% =========== Interpolate velocity onto the regular grid =========== %
    [ID]        =   InterpStaggered(D,ID,N,'velocity');
    D.meanV(it) =   rms(ID.vx(:) + ID.vz(:));
    % =================================================================== %
    %% ========================== Plot data ============================= %
    Pl.time     =   ...
        ['@ Iteration: ',sprintf('%i',it),...
        '; Time: ',sprintf('%2.2e',T.time(it)/1e6/(365.25*24*60*60)),' Myr'];
    
    switch Pl.plotfields
        case 'yes'
            if (mod(it,5)==0||it==1)
                figure(3)
                clf
                
                ax1=subplot(3,1,1);
                plotfield(D.rho,M.X/1e3,M.Z/1e3,Pl,'contourf',...
                    '\it\rho\rm\bf','quiver',ID.vx,ID.vz)
                colormap(ax1,flipud(Pl.oslo))
                ax2=subplot(3,1,2);
                plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
                    '\it\eta\rm\bf','contour2',D.C,1.5)
                colormap(ax2,flipud(Pl.lapaz))
                ax3=subplot(3,1,3);
                plotfield(ID.v.*100*365.25*24*60*60,...
                    M.X/1e3,M.Z/1e3,Pl,'pcolor','\itv \rm\bf')
                colormap(ax3,Pl.imola)
                
                switch Pl.savefig
                    case 'yes'
                        saveas(figure(3),...
                            [M.ModDir,'/Field',num2str(it)],'png')
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
            end
        case 'no'
            disp(['Iteration: ',sprintf('%i',it)])
    end
    % =================================================================== %
    %% ========================== Advection ============================= %
    [D,Ma,ID]       =   Advection(it,N,B,D,ID,Py,T.dt,M,Ma);
    % =================================================================== %
end
%% ======================== Save final figure =========================== %
figure(3)
clf
ax1=subplot(3,1,1);
plotfield(D.rho,M.X/1e3,M.Z/1e3,Pl,'contourf',...
    '\it\rho\rm\bf','quiver',ID.vx,ID.vz)
colormap(ax1,flipud(Pl.oslo))
ax2=subplot(3,1,2);
plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
    '\it\eta\rm\bf')
colormap(ax2,flipud(Pl.lapaz))
ax3=subplot(3,1,3);
plotfield(ID.v.*100*365.25*24*60*60,...
    M.X/1e3,M.Z/1e3,Pl,'pcolor','\itv \rm\bf')
colormap(ax3,Pl.imola)
switch Pl.savefig
    case 'yes'
        saveas(figure(3),['data/RTI_adv_',B.AdvMethod,...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nz_',num2str(N.nz),'_SS'],'png')
end
% ======================================================================= %
%% ======================== Plot time serieses ========================== %
figure(4)
plot(T.time/1e6/(365.25*24*60*60),...
    D.meanV*100*365.25*24*60*60,...
    'LineWidth',2)
set(gca,'FontWeight','Bold');
xlabel('t [ Ma ]'); ylabel('VRMS [ cm/a ]')
title('Root Mean Square Velocity')

switch Pl.savefig
    case 'yes'
        saveas(figure(4),[M.ModDir,'/TimeSeries'],'png')
end
% ======================================================================= %
T.tend      = toc(T.tstart);
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


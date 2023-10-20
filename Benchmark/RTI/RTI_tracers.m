clear
clc
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\AdvectionProblem')
    addpath('..\..\StokesProblem')
    addpath('..\..\SetUp')
else
    addpath('../../AdvectionProblem')
    addpath('../../StokesProblem')
    addpath('../../SetUp')
end
% ======================================================================= %
%% ===================== Some initial definitions ======================= %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'tracers';
B.Aparam        =   'comp';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'variable';
B.EtaIni        =   'RTI';
B.lambda        =   0.5; 		% Wellenlaenge der Perturbation [ km ]
B.deltaA        =   100;        % Amplitude [ km ]
% ======================================================================= %
%% ========================= Define flow field ========================== %
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
T.itmax     =   50;          %   Maximal erlaubte Anzahl der Iterationen
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
Pl.xlab     =   'x [km]';
Pl.zlab     =   'z [km]';
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
%% ================ Information for the command window ================== %
fprintf([' Rayleigh-Taylor Instabilitaet  --------------------- ',...
    '\n Advektion mit: %s',...
    '\n Viskositaet ist: %s',...
    '\n Aufloesung (nx x nz): %i x %i',...
    '\n Refrenzviskositaet [Pa s]: %2.2e',...
    '\n  --------------------- ',...
    '\n\n '],B.AdvMethod,Py.eparam,...
    N.nx,N.nz,Py.eta0);
fprintf(['Maximum Time : %1.4g',...
    '\n Maximale Anzahl an Iterationen: %5i',...
    '\n  --------------------- \n'],...
    T.tmax/1e6/(365.25*24*60*60),T.itmax)
% ======================================================================= %
%% ========================= Time loop ================================= %%
for it = 1:T.itmax
    if(strcmp(B.AdvMethod,'none')==0)
        [D,A]       =   solveSECE(D,Py,N,B);
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
        ['Time: ',sprintf('%2.2e',T.time(it)/1e6/(365.25*24*60*60)),...
        ' Myr'];
    
    switch Pl.plotfields
        case 'yes'
            if (mod(it,5)==0||it==1)
                figure(3)
                clf
                ax1=subplot(3,1,1);
                Pl.cbtitle  =   [{'$$\rho$$'},{'$$[kg/m^3]$$'}];
                plotfield(D.rho,M.X/1e3,M.Z/1e3,Pl,'contourf',...
                    [],'quiver',ID.vx,ID.vz)
                colormap(ax1,flipud(Pl.oslo))
                ax2=subplot(3,1,2);
                Pl.cbtitle  =   [{'$$\eta$$'},{'$$[ Pa s ]$$'}];
                plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
                    [])
                colormap(ax2,flipud(Pl.lapaz))
                ax3=subplot(3,1,3);
                Pl.cbtitle  =   [{'v'},{'[ cm/a ]'}];
                plotfield(ID.v.*100.*(60*60*24*365.25),...
                    M.X/1e3,M.Z/1e3,Pl,'pcolor',[])
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
Pl.cbtitle  =   [{'$$\rho$$'},{'$$[kg/m^3]$$'}];
plotfield(D.rho,M.X/1e3,M.Z/1e3,Pl,'contourf',...
    [],'quiver',ID.vx,ID.vz)
colormap(ax1,flipud(Pl.oslo))
ax2=subplot(3,1,2);
Pl.cbtitle  =   [{'$$\eta$$'},{'$$[ Pa s ]$$'}];
plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
    [])
colormap(ax2,flipud(Pl.lapaz))
ax3=subplot(3,1,3);
Pl.cbtitle  =   [{'v'},{'[ cm/a ]'}];
plotfield(ID.v.*100.*(60*60*24*365.25),...
    M.X/1e3,M.Z/1e3,Pl,'pcolor',[])
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
set(gca,'FontWeight','Bold',...
    'LineWidth',2,'FontSize',15,'TickLabelInterpreter','latex')
xlabel('$$t\ [Ma]$$','Interpreter','latex')
ylabel('$$V_{RMS}\ [cm/a]$$','Interpreter','latex')
title('$$Root\ Mean\ Square\ Velocity$$','Interpreter','latex')

switch Pl.savefig
    case 'yes'
        saveas(figure(4),[M.ModDir,'/TimeSeries'],'png')
end
% ======================================================================= %
%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\..\AdvectionProblem')
    rmpath('..\..\StokesProblem')
    rmpath('..\..\SetUp')
else
    rmpath('../../AdvectionProblem')
    rmpath('../../StokesProblem')
    rmpath('../../SetUp')
end
% ======================================================================= %
% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


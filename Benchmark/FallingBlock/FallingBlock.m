clear
clc
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
% ======================================================================= %
T.tstart        =   tic;
%% ===================== Some initial definitions ======================= %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
Py.scale        =   'no';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'tracers';
B.Aparam        =   'comp';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'variable';
B.EtaIni        =   'block';
% ======================================================================= %
%% ========================= Define flow field ========================== %
B.FlowFac       =   10;
% ======================================================================= %
%% ==================== Define model geometry constants ================= %
M.H         =   -500;          %   Depth [ in km ]
M.xmax      =   1;              %   Aspect ratio
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   51;             %   Vertical grid solution
N.nx        =   51;             %   Horizontal grid solution
% ======================================================================= %
%% ====================== Tracer advection method ======================= %
% Number of tracers per cell: nmx*nmz                                     %
N.nmx       =   5;
N.nmz       =   5;
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   9.81;                 %   Gravitational acceleration [m/s^2]
Py.rho0     =   3200;               %   Reference density [kg/m^3]
Py.k        =   5;                  %   Thermal conductivity [ W/m/K ]
Py.cp       =   1250;               %   Heat capacity [ J/kg/K ]
Py.alpha    =   2.5e-5;             %   Thermal expansion coeff. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermal diffusivity [ m^2/s ]

Py.Q0       =   0; %3.1e-9;         %   Heat production rate per volume [W/m^3]
Py.Q0       =   Py.Q0/Py.rho0;      %   Heat production rate per mass [W/kg]

Py.eta0     =   1e21;               %   Viscosity comp. 0 [ Pa*s ]
Py.eta1     =   1e21;               %   Viscosity comp. 1 [ Pa*s ]

Py.rho1     =   3300;               %   Density comp. 1 [ kg/m^3 ]

Py.DeltaT   =   1000;               %   Temperature difference top-bottom
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Velocity boundary conditions ------------------------------------------ %
%   0 - no slip; 1 - free slip;
B.tbc       =   1;              %   Top
B.bbc       =   1;              %   Bottom
B.lbc       =   1;              %   Left
B.rbc       =   1;              %   Right
% ======================================================================= %
%% ====================== Define time constants ========================= %
T.tmaxini   =   9.886;             %   Maximum time in Ma
T.itmax     =   50;            %   Maximum number of iterations
T.dtfac     =   0.5;            %   Courant criterium
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
Pl.xlab     =   '\bfx';
Pl.zlab     =   '\bfz';
switch Pl.plotfields
    case 'yes'
        if strcmp(getenv('OS'),'Windows_NT')
            set(figure(2),'position',[1.8,1.8,766.4,780.8]);
            h           =   figure(2);
        else
            set(figure(2),'position',[-1919,1,960,988]);
            h           =   figure(2);
        end
end
% Animation settings ---------------------------------------------------- %
switch Pl.savefig
    case 'yes'
        M.ModDir    = ['data/Blanckenbach_Ra_',sprintf('%2.2e',Py.Ra),...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nx_',num2str(N.nz),'_nz_',num2str(N.nz)];
        if ~exist(M.ModDir,'dir')
            mkdir(M.ModDir)
        end
        Pl.filename    = [M.ModDir,'/Evolution.gif'];
        set(figure(2),'position',[1.8,1.8,766.4,780.8]);
        h           =   figure(2);
end
% ======================================================================= %
%% ========================== Scale Parameters ========================== %
switch lower(Py.scale)
    case 'yes'
        [M,N,D,T,S]     =   ScaleParameters(M,Py,N,D,T);
end
% ======================================================================= %
%% ================ Information for the command window ================== %
fprintf([' Falling Block --------------------- ',...
    '\n Advektion mit: %s',...
    '\n Viskositaet ist: %s',...
    '\n Aufloesung (nx x nz): %i x %i',...
    '\n Refrenzviskositaet [Pa s]: %2.2e',...
    '\n  --------------------- ',...
    '\n '],B.AdvMethod,Py.eparam,...
    N.nx,N.nz,Py.eta0);
fprintf(['Maximum Time : %1.4g',...
    '\n Maximale Anzahl an Iterationen: %5i',...
    '\n  --------------------- \n'],...
    T.tmax,T.itmax)
% ======================================================================= %
%% ========================= Time loop ================================= %%
for it = 1:T.itmax
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
    if T.time(it) > T.tmax
        T.dt        =   T.tmax - T.time(it-1); 
        T.time(it)  =   T.time(it-1) + T.dt; 
    end
    % =================================================================== %
    %% =========== Interpolate velocity onto the regular grid =========== %
    [ID]        =   InterpStaggered(D,ID,N,'velocity');
    D.meanV(it) =   mean(ID.v,'all');   % Mittleregeschwindigkeit
    D.meanV(it) =   rms(ID.vx(:) + ID.vz(:));
    % =================================================================== %
    %% ========================== Plot data ============================= %
    Pl.time     =   ...
        ['@ Iteration: ',sprintf('%i',it),...
        '; Time: ',sprintf('%2.2e',T.time(it))];
    if (mod(it,1)==0||it==1)
        switch Pl.plotfields
            case 'yes'
                figure(2)
                clf
                switch Py.eparam
                    case 'const'
                        ax1=subplot(2,1,1);
                        plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                            '\it\rho \rm\bf','quiver',ID.vx,ID.vz)
                        colormap(ax1,flipud(Pl.lajolla))
                        ax2=subplot(2,1,2);
                        plotfield(ID.v,M.X./1e3,M.Z/1e3,Pl,'pcolor',...
                            '\itv \rm\bf')
                        colormap(ax2,Pl.imola)
                    case 'variable'
                        ax1=subplot(2,1,1);
                        plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                            '\it\rho\rm\bf','quiver',ID.vx,ID.vz)
                        colormap(ax1,flipud(Pl.lajolla))
                        ax2=subplot(2,1,2);
                        plotfield(log10(D.eta),M.X./1e3,M.Z./1e3,Pl,'contourf',...
                            '\it\eta\rm\bf','quiver',ID.vx,ID.vz)
                        colormap(ax2,flipud(Pl.lapaz))
                end
                switch Pl.savefig
                    case 'yes'
                        saveas(figure(2),...
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
            case 'no'
                disp(['Iteration: ',sprintf('%i',it)])
        end
    end
    % =================================================================== %
    %% ========================== Advection ============================= %
    [D,Ma,ID]       =   Advection(it,B,D,ID,Py,T.dt,M,Ma);
    % =================================================================== %
    %% ========================== Check break =========================== %
    if (T.time(it) >= T.tmax)
        set(figure(3),'position',[1.8,1.8,766.4,780.8]);
        h           =   figure(3);
        disp('Maximum time reached, stop time loop!')
        T.indtime   = find(T.time(2:end)==0,1,'first');
        figure(3)
        clf
        switch Py.eparam
            case 'const'
                ax1=subplot(2,1,1);
                plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                    '\it\rho \rm\bf','quiver',ID.vx,ID.vz)
                colormap(ax1,flipud(Pl.lajolla))
                ax2=subplot(2,1,2);
                plotfield(ID.v,M.X./1e3,M.Z/1e3,Pl,'pcolor',...
                    '\itv \rm\bf')
                colormap(ax2,Pl.imola)
            case 'variable'
                ax1=subplot(2,1,1);
                plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                    '\it\rho\rm\bf','quiver',ID.vx,ID.vz)
                colormap(ax1,flipud(Pl.lajolla))
                ax2=subplot(2,1,2);
                plotfield(log10(D.eta),M.X./1e3,M.Z./1e3,Pl,'contourf',...
                    '\it\eta\rm\bf','quiver',ID.vx,ID.vz)
                colormap(ax2,flipud(Pl.lapaz))
        end
        break
    end
    % =================================================================== %
end
%% ======================== Save final figure =========================== %
switch Pl.savefig
    case 'yes'
        saveas(figure(3),[M.ModDir,'/Field_SS'],'png')
end
% ======================================================================= %
T.tend      = toc(T.tstart);
%% ====================== Clear path structure ========================== %
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
% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


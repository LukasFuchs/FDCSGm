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
Py.scale        =   'no';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'variable';
B.EtaIni        =   'block';
% ======================================================================= %
%% ========================= Define flow field ========================== %
B.FlowFac       =   10;
% ======================================================================= %
eta0            =   1e21;
eta1            =   logspace(15,27,13);
vmax            =   zeros(1,length(eta1));
etar            =   zeros(1,length(eta1));
for k = 1:length(eta1)
    %% ==================== Define model geometry constants ============= %
    M.H         =   -500;          %   Depth [ in km ]
    M.xmax      =   1;              %   Aspect ratio
    % =================================================================== %
    %% ====================== Define the numerical grid ================= %
    N.nz        =   51;             %   Vertical grid solution
    N.nx        =   51;             %   Horizontal grid solution
    % =================================================================== %
    %% ====================== Tracer advection method =================== %
    % Number of tracers per cell: nmx*nmz                                 %
    N.nmx       =   5;
    N.nmz       =   5;
    % =================================================================== %
    %% ====================== Define physical constants ================= %
    Py.g        =   9.81;               %   Gravitational acceleration [m/s^2]
    Py.rho0     =   3200;               %   Reference density [kg/m^3]
    Py.k        =   5;                  %   Thermal conductivity [ W/m/K ]
    Py.cp       =   1250;               %   Heat capacity [ J/kg/K ]
    Py.alpha    =   2.5e-5;             %   Thermal expansion coeff. [ K^-1 ]
    
    Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermal diffusivity [ m^2/s ]
    
    Py.Q0       =   0; %3.1e-9;         %   Heat production rate per volume [W/m^3]
    Py.Q0       =   Py.Q0/Py.rho0;      %   Heat production rate per mass [W/kg]
    
    Py.eta0     =   eta0;               %   Viscosity comp. 0 [ Pa*s ]
    Py.eta1     =   eta1(k);            %   Viscosity comp. 1 [ Pa*s ]
    
    Py.rho1     =   3300;               %   Density comp. 1 [ kg/m^3 ]
    
    Py.DeltaT   =   1000;               %   Temperature difference top-bottom
    % =================================================================== %
    %% ===================== Define boundary conditions ================= %
    % Velocity boundary conditions -------------------------------------- %
    %   0 - no slip; 1 - free slip;
    B.tbc       =   1;              %   Top
    B.bbc       =   1;              %   Bottom
    B.lbc       =   1;              %   Left
    B.rbc       =   1;              %   Right
    % =================================================================== %
    %% ====================== Define time constants ===================== %
    T           =   [];
    % =================================================================== %
    %% ========================= Define fields required ================= %
    [Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
    % =================================================================== %
    %% ======================== Setup initial conditions ================ %
    [T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
    % =================================================================== %
    %% ========================= Plot parameter ========================= %
    Pl.inc      =   min(N.nz/10,N.nx/5);
    Pl.inc      =   round(Pl.inc);
    Pl.xlab     =   '$$x$$';
    Pl.zlab     =   '$$z$$';
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
    switch Pl.savefig
        case 'yes'
            M.ModDir    = ['data/FallingBlock',...
            '_etar_',num2str(Py.eta1/Py.eta0),...
            '_drho_',num2str(Py.rho1-Py.rho0),...
            '_nx_',num2str(N.nz),'_nz_',num2str(N.nz)];
            if ~exist(M.ModDir,'dir')
                mkdir(M.ModDir)
            end
    end
    % =================================================================== %
    %% ================ Information for the command window ============== %
    fprintf([' Falling Block --------------------- ',...
        '\n Viskositaet ist: %s',...
        '\n Aufloesung (nx x nz): %i x %i',...
        '\n Refrenzviskositaet [Pa s]: %2.2e',...
        '\n  --------------------- ',...
        '\n '],Py.eparam,...
        N.nx,N.nz,Py.eta0);
    % =================================================================== %
    %% ============================ Solving equations =================== %
    switch Py.eparam
        case 'const'
            [D,A]       =   solveSECE_const_Eta(D,Py,N,B,A);
        case 'variable'
            [D,A]       =   solveSECE(D,Py,N,B);
        otherwise
            error('Viscosity not difined! Check Py.eparam parameter.')
    end
    % =================================================================== %
    %% =========== Interpolate velocity onto the regular grid =========== %
    [ID]        =   InterpStaggered(D,ID,N,'velocity');
    % =================================================================== %
    %% ========================== Plot data ============================= %
    Pl.time     =   '';
    switch Pl.plotfields
        case 'yes'
            figure(2)
            clf
            switch Py.eparam
                case 'const'
                    ax1=subplot(2,1,1);
                    plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
                        '$$\rho$$','quiver',ID.vx,ID.vz)
                    colormap(ax1,flipud(Pl.lajolla))
                    ax2=subplot(2,1,2);
                    plotfield(ID.v,M.X./1e3,M.Z/1e3,Pl,'pcolor',...
                        '$$v$$')
                    colormap(ax2,Pl.imola)
                case 'variable'
                    ax1=subplot(2,1,1);
                    plotfield(D.rho,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                        '$$\rho$$','quiver',ID.vx,ID.vz)
                    colormap(ax1,flipud(Pl.lajolla))
                    ax2=subplot(2,1,2);
                    plotfield(log10(D.eta),M.X./1e3,M.Z./1e3,Pl,'pcolor',...
                        '$$\eta$$')
                    colormap(ax2,flipud(Pl.lapaz))
            end
            switch Pl.savefig
                case 'yes'
                    saveas(figure(2),...
                        [M.ModDir,'/FieldSS'],'png')
            end
        case 'no'
            disp(['Iteration: ',sprintf('%i',it)])
    end
    % =================================================================== %
    vmax(k)     =   max(ID.v(M.ind));
    etar(k)     =   (Py.eta1/Py.eta0);
end
figure(3)
plot(etar,vmax,'kd-','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k')
xlabel('$$log_{10}( \eta_{block} / \eta_{medium} )$$','Interpreter','latex')
ylabel('Block velocity [ m/s ]','Interpreter','latex')
set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,'xscale','log',...
    'TickLabelInterpreter','latex')
switch Pl.savefig
    case 'yes'
        saveas(figure(3),['data',...
                '/drho_',num2str(Py.rho1-Py.rho0),...
                '_vz_eta_r_nx_',num2str(N.nx)],'png')
end
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


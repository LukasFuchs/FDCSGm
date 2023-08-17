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
T.tstart        =   tic;
%% ===================== Some initial definitions ======================= %
% Save figures - 'yes' or 'no'                                            %
Pl.savefig      =   'no';

% Plot fields - 'yes' or 'no'                                             %
Pl.plotfields   =   'yes';

% Scale parameters - 'yes' or 'no'                                        %
Py.scale        =   'no';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
% Choose advection method:                                                %
%   'semi-lag'                                                            %
%   'Upwind'                                                              %
%   'SLF'                                                                 %
%   'tracers'                                                             %
B.AdvMethod     =   'semi-lag';

% Define which property should be advected:                               %
%   'none' 								  %
%   'rho' 								  %
%   'comp'                                                                %
%   'temp' 								  %
B.Aparam        =   'comp';

% Choose T-diffusion method:                                              %
%   'none'                                                                %
%   'ADI'                                                                 %
%   'CNV'                                                                 %
%   'explicit'                                                            %
%   'implicit'                                                            %
B.DiffMethod    =   'ADI';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
% Switch to define viscosity field:                                       %
%   'none'        - no stokes solution                                    %
%   'const'                                                               %
%   'variable'                                                            %
Py.eparam       =   'const';

% Define initial viscosity field ---------------                          %
%   'none' 								  %
%   'tdep'                                                                %
%   'ellipse'                                                             %
%   'RTI'                                                                 %
B.EtaIni        =   'tdep';

% If B.EtaIni == 'tdep' ------------------------                          %
% Define parameters for T-dep viscosity                                   %
%   1) b, c -> eta = eta0 * exp(-b*T/dT + c*z/H)                          %
Py.b            =   log(1000); %log(16384);     % Temperaturabhaengigkeit
Py.c            =   0; %log(64);                % Tiefenabhaengigkeit

% If EtaIni == 'ellipse' ----------------------                           %
%   Define parameters of the ellipse:                                     %
%       RotAng  -   Rotation angle; positive -> counter clockwise         %
%       EllA    -   Major-half axis in [ m ]                              %
%       EllB    -   Minor-half axis in [ m ]                              %
B.RotAng        =   0;
B.EllA          =   2e2;
B.EllB          =   2e2;

% If EtaIni == 'RTI' ------------------------                             %
%   Define parameters of the initial perturbation:                        %
%       lambda  -   Wavelength in km                                      %
%       deltaA  -   Amplitude in km                                       %
B.lambda        =   0.5;
B.deltaA        =   100;
% ======================================================================= %
%% ================== Define initial temperature anomaly ================ %
%   'none' 								  %
%   'ellipse'                                                             %
%   'gaussian'                                                            %
%   'block'                                                               %
%   'const'                                                               %
%   'linear'                                                              %
%   'linano'                                                              %
B.Tini          =   'block';

% If Tini == 'gaussian' ----------------------                            %
%   Define parameters for the gaussian temperature distribution:          %
%       Ampl    -   Amplitude [ K ]                                       %
%       sigma   -   Width; percentage of L [ m ]                          %
%       T0      -   Background temperature [ K ]                          %
%   The gaussian is by default put at a position of 1/4*L and 0.5*H.      %
B.TAmpl         =   10;
B.Tsigma        =   10;
B.T0            =   1000;

% Define if thermal parameters are 'const' or 'variable' or 'none'        %
%   'none' 	- No energy equation solution 				  %
%   'const' 	- 						 	  %
%   'variable'  - 							  %
Py.tparam       =   'const';
% ======================================================================= %
%% ========================= Define flow field ========================== %
%   'none'                                                                %
%   'RigidBody'                                                           %
%   'ShearCell'                                                           %
%   'PureShear'                                                           %
%   'SimpleShear'                                                         %
B.IniFlow       =   'none';

% If B.IniFlow == 'PureShear' || 'SimpleShear'                            %
%   Define background strain rate ebg [ s^-1 ]; % < 0 compression         %
B.ebg           =   -1e15;

B.FlowFac       =   10;
% ----------------------------------------------------------------------- %
%% ==================== Define model geometry constants ================= %
M.H         =   -2900;          %   Depth [ in km ]
M.xmax      =   1;              %   Aspect ratio

M.zUC       =   -10e3;          %   Depth of the upper crust [ km ]
M.zLC       =   -35e3;          %   Depth of the lower crust [ km ]
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
Py.g        =   10;                 %   Gravitational acceleration [m/s^2]
Py.rho0     =   4000;               %   Reference density [kg/m^3]
Py.k        =   5;                  %   Thermal conductivity [ W/m/K ]
Py.cp       =   1250;               %   Heat capacity [ J/kg/K ]
Py.alpha    =   2.5e-5;             %   Thermal expansion coeff. [ K^-1 ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermal diffusivity [ m^2/s ]

Py.Q0       =   7.45e-09; %3.1e-9;  %   Heat production rate per volume [W/m^3]
Py.Q0       =   Py.Q0/Py.rho0;      %   Heat production rate per mass [W/kg]

Py.eta0     =   1e23;               %   Viscosity comp. 0 [ Pa*s ]
Py.eta1     =   1e23;               %   Viscosity comp. 1 [ Pa*s ]

Py.rho1     =   3000;               %   Density comp. 1 [ kg/m^3 ]

Py.DeltaT   =   1000;               %   Temperature difference top-bottom

% If Ra < 0, then Ra is calculated from the parameters above. If Ra is    %
% given, then the reference viscosity eta0 is adjusted to obtain Ra using %
% the other given scaling parameters. Ra is defined as:                   %
%   Ra  = (rho0 * g * alpha * DeltaT * H^3) / (eta0*kappa)                %
Py.Ra       =   1e6;                %   Rayleigh number
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Velocity boundary conditions ------------------------------------------ %
%   0 - no slip; 1 - free slip;
B.tbc       =   1;              %   Top
B.bbc       =   1;              %   Bottom
B.lbc       =   1;              %   Left
B.rbc       =   1;              %   Right

% Thermal boundary conditions ------------------------------------------- %
%   'const' -
%   'flux'  -
B.ttbc      =   'const';        %   Top
B.btbc      =   'flux';         %   Bottom
B.ltbc      =   'flux';         %   Left
B.rtbc      =   'flux';         %   Right
% Define actual temperature condition:
%   if 'flux'   -   following parameters define the flux through the boundary
%   if 'const'  -   following parameters define the temperatue at the boundary in K
B.lhf       =   0;
B.rhf       =   0;
B.thf       =   273;
B.bhf       =   B.thf + Py.DeltaT;
% ======================================================================= %
%% ====================== Define time constants ========================= %
T.tmaxini   =   17000;          %   Maximum time in Ma
T.itmax     =   5000;           %   Maximum number of iterations
T.dtfac     =   1.0;            %   Courant criterium
T.dtdifac   =   0.8;            %   Diffusion stability criterium
% ======================================================================= %
%% ========================= Define fields required ===================== %
[Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
% ======================================================================= %
%% ======================== Setup initial conditions ==================== %
[T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
% ======================================================================= %
%% ======================= Rayleigh number conditions =================== %
if Py.Ra < 0
    % Falls die Rayleigh Zahl nicht explizit angegeben wird, wird sie hier
    % berechnet.
    Py.Ra   =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.eta0/Py.kappa;
else
    % Falls die Rayleigh Zahl gegeben ist m??ssen wir eine Variable
    % anpassen, z.B. die Referenzviskosit??t.
    Py.eta0 =   Py.rho0*Py.g*Py.alpha*Py.DeltaT*(-M.H)^3/Py.Ra/Py.kappa;
end
% ======================================================================= %
%% ========================= Plot parameter ============================= %
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
        M.ModDir    = ['data/Blanckenbach_Ra_',sprintf('%2.2e',Py.Ra),...
            '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
            '_nx_',num2str(N.nz),'_nz_',num2str(N.nz)];
        if ~exist(M.ModDir,'dir')
            mkdir(M.ModDir)
        end
        Pl.filename    = [M.ModDir,'/Evolution.gif'];
        set(figure(1),'position',[1.8,1.8,766.4,780.8]);
        h           =   figure(1);
end
% ======================================================================= %
%% ========================== Scale Parameters ========================== %
switch lower(Py.scale)
    case 'yes'
        [M,N,D,T,S]     =   ScaleParameters(M,Py,N,D,T);
end
% ======================================================================= %
%% ================ Information for the command window ================== %
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
% ======================================================================= %
%% ========================= Time loop ================================= %%
for it = 1:T.itmax
    %     disp(['Iteration: ',sprintf('%i',it)])
    if(strcmp(B.AdvMethod,'none')==0)
        switch Py.eparam
            case 'none'
                disp('No stokes solution activated!')
            case 'const'
                switch Py.scale
                    case 'yes'
                        [D,A]       =   solveSECE_const_EtaSc(D,Py,N,B,A);
                    case 'no'
                        [D,A]       =   solveSECE_const_Eta(D,Py,N,B,A);
                end
                if (it == 2)
                    N.beenhere = 1;
                end
            case 'variable'
                switch Py.scale
                    case 'yes'
                        [D,A]       =   solveSECESc(D,Py,N,B);
                    case 'no'
                        [D,A]       =   solveSECE(D,Py,N,B);
                end
            otherwise
                error('Viscosity not difined! Check Py.eparam parameter.')
        end
    end
    % =================================================================== %
    %% ===================== Calculate time stepping ==================== %
    T.dt        =   T.dtfac*min(N.dx,abs(N.dz))/...
        (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));
    switch Py.scale
        case 'no'
            T.dtdiff    =   T.dtdifac*min([N.dx,abs(N.dz)])^2/Py.kappa/4;
        case 'yes'
            T.dtdiff    =   T.dtdifac*min([N.dx,abs(N.dz)])^2/4;
    end
    T.dt        =   min(T.dt,T.dtdiff);
    if it>1
        T.time(it)  =   T.time(it-1) + T.dt;
    end
    % =================================================================== %
    %% =========== Interpolate velocity onto the regular grid =========== %
    [ID]        =   InterpStaggered(D,ID,N,'velocity');
    D.meanV(it) =   mean(ID.v,'all');   % Mittleregeschwindigkeit
    D.meanV(it) =   rms(ID.vx(:) + ID.vz(:));
    % =================================================================== %
    %% =================== Calculate strain rate ======================== %
    [ID]        =   GetStrainRate(ID,N);
    % =================================================================== %
    %% ========================== Plot data ============================= %
    Pl              =   PlotData(it,Pl,T,D,M,ID,Py);
    % Alternatively like below; similar plotting routines should be       %
    % included in the PlotData function from above! 			  %
    Pl.time     =   ...
        ['@ Iteration: ',sprintf('%i',it),...
        '; Time: ',sprintf('%2.2e',T.time(it))];
    if (mod(it,10)==0||it==1)
        switch Pl.plotfields
            case 'yes'
                figure(1)
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
    % Alternatively like below; the same routines should be included in   %
    % the Advection function above! 				  	  %
    switch lower(B.AdvMethod)
        case 'tracers'
            % Advect tracers with Runge-Kutta 4th order ------------- %
            [Ma.XM,Ma.ZM] = ...
                AdvectMarker2D(M,Ma.XM,Ma.ZM,dt,ID.vx,ID.vz);
            % Interpolate from the tracers to the grid -------------- %
            % This is so far only recommendable if eta and rho are    %
            % constant for each composition. Needs some modification  %
            % to make eta tdep etc.					  %
            [~,D.rho]   =   TracerInterp(Ma,D.rho,[],M.X,M.Z,'from','rho');
            [~,D.eta]   =   TracerInterp(Ma,D.eta,[],M.X,M.Z,'from','eta');
        case 'semi-lag'
            switch lower(B.Aparam)
                case 'rho'
                    if (it==1)
                        [D.rho]     =   SemiLagAdvection2D(ID,M,D.rho,dt);
                        
                        % Speicher die alte Geschwindigkeit
                        ID.vxo      =   ID.vx;
                        ID.vzo      =   ID.vz;
                    else
                        [D.rho]       =   SemiLagAdvection2D(ID,M,D.rho,dt);
                    end
                case 'temp'
                    if (it==1)
                        [D.T]     =   SemiLagAdvection2D(ID,M,D.T,dt);
                        
                        % Speicher die alte Geschwindigkeit
                        ID.vxo          =   ID.vx;
                        ID.vzo          =   ID.vz;
                    else
                        [D.T]     =   SemiLagAdvection2D(ID,M,D.T,dt);
                    end
                case 'comp'
                    if (it==1)
                        [D.C]     =   SemiLagAdvection2D(ID,M,D.C,dt);
                        
                        % Speicher die alte Geschwindigkeit
                        ID.vxo          =   ID.vx;
                        ID.vzo          =   ID.vz;
                    else
                        [D.C]     =   SemiLagAdvection2D(ID,M,D.C,dt);
                    end
                    % Update viscosity and density ========================== %
                    if (mod(it,25)==0)
                        D.C             =   round(D.C);
                    end
                    D.rho(D.C<=1.5) =    Py.rho0;   D.rho(D.C>=1.5) =    Py.rho1;
                    D.eta(D.C<=1.5) =    Py.eta0;   D.eta(D.C>=1.5) =    Py.eta1;
            end
    end
end
% =================================================================== %
%% ========================== Diffusion ============================= %
if ~strcmp(B.Tini,'none')
    D.T             =   Diffusion(B,D.T,D.Q,D.rho,Py,T.dt,N);
end
% Alternatively like this; the same routines should be included in    %
% the Diffusion function above.                                       %
if isfield(B,'DiffMethod')
    switch Py.scale
        case 'no'
            switch lower(B.DiffMethod)
                case 'explicit'
                    D.T     =   SolveDiff2Dexplicit(D.T,D.Q,D.rho,T.dt,Py,N,B);
                case 'implicit'
                    D.T     =   SolveDiff2Dimplicit(D.T,D.Q,D.rho,T.dt,Py,N,B);
                case 'adi'
                    D.T     =   SolveDiff2DADI(D.T,D.Q,D.rho,T.dt,Py,N,B);
                case 'cnv'
                    D.T     =   SolveDiff2DCNV(D.T,D.Q,D.rho,T.dt,Py,N,B);
                case 'none'
                otherwise
                    error('Diffusion scheme is not defined!')
            end
        case 'yes'
            switch lower(B.DiffMethod)
                case 'explicit'
                    D.T     =   SolveDiff2DexplicitSc(D.T,D.Q,T.dt,N,B);
                case 'implicit'
                    D.T     =   SolveDiff2DimplicitSc(D.T,D.Q,T.dt,N,B);
                    %         case 'implicito'
                    %             [D.T,N] =   SolveDiff2Dimplicit_optSc(D.T,D.Q,T.dt,N,B);
                case 'adi'
                    D.T     =   SolveDiff2DADISc(D.T,D.Q,T.dt,N,B);
                case 'none'
                otherwise
                    error('Diffusion scheme is not defined!')
            end
    end
end
% =================================================================== %
%% ================== Heat flow at the surface ====================== %
D.dTtop         =   (D.T(1,:)-D.T(2,:))./N.dz;
D.dTbot         =   (D.T(end-1,:)-D.T(end,:))./N.dz;
%     D.Nus(it)       =   trapz(M.x,D.dTtop).*abs(M.H)/M.L;
D.Nus(it)       =   mean(D.dTtop);

D.meanT(:,it)   =   mean(D.T,2);
% =================================================================== %
%% ====================== Update viscosity ========================== %
switch lower(B.EtaIni)
    case 'tdep'
        D.eta   =   exp( -Py.b.*((D.T-D.T(1,1))./(D.T(N.nz,1)-D.T(1,1)))...
            + Py.c.*M.Z);
end
% =================================================================== %
%% ========================== Check break =========================== %
[answer,T]  =   CheckBreakCriteria(it,T,D,M,Pl,ID,Py);
switch answer
    case 'yes'
        break
end
% =================================================================== %
end
%% ======================== Save final figure =========================== %
switch Pl.savefig
    case 'yes'
        saveas(figure(1),[M.ModDir,'/Field_SS'],'png')
end
% ======================================================================= %
%% ======================== Plot time serieses ========================== %
PlotTimeSerieses(Py,T,D,M,N)
switch Pl.savefig
    case 'yes'
        saveas(figure(2),[M.ModDir,'/TimeSeries'],'png')
end
% ======================================================================= %
T.tend      = toc(T.tstart);
%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\DiffusionProblem')
    rmpath('..\AdvectionProblem')
    rmpath('..\StokesProblem')
    rmpath('..\SetUp')
    rmpath('..\ScaleParam')
else
    rmpath('../DiffusionProblem')
    rmpath('../AdvectionProblem')
    rmpath('../StokesProblem')
    rmpath('../SetUp')
    rmpath('../ScaleParam')
end
% ======================================================================= %
profile viewer
% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


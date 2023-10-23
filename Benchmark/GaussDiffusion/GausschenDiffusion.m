clear
clc
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\DiffusionProblem')
    addpath('..\..\SetUp')
    addpath('..\..\ScaleParam')
else
    addpath('../../DiffusionProblem')
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
B.DiffMethod    =   'explicit';
% ======================================================================= %
%% ================== Define initial temperature anomaly ================ %
B.Tini          =   'gaussian';
B.TAmpl         =   200;            % Temperature amplitude [ K ]
B.Tsigma        =   10;             % Percentage of L [ % ]
B.T0            =   1000;           % Background temperature [ K }

Py.tparam       =   'const';
% ======================================================================= %
%% ==================== Define model geometry constants ================= %
M.H         =   -200;               %   Depth [ in km ]
M.xmax      =   1;                  %   Aspect ratio
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   51;                 %   Vertical grid solution
N.nx        =   51;                 %   Horizontal grid solution
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   10;                 %   Gravitational acceleration [m/s^2]
Py.rho0     =   3200;               %   Reference density [kg/m^3]
Py.k        =   3;                  %   Thermal conductivity [ W/m/K ]
Py.cp       =   1000;               %   Heat capacity [ J/kg/K ]

Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermal diffusivity [ m^2/s ]
Py.DeltaT   =   B.T0;
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
B.ttbc      =   'const';         %   Top
B.btbc      =   'const';         %   Bottom
B.ltbc      =   'const';         %   Left
B.rtbc      =   'const';         %   Right
% Define actual temperature condition:
%   if 'flux'   -   following parameters define the flux through the boundary
%   if 'const'  -   following parameters define the temperatue at the boundary in K
B.lhf       =   B.T0;
B.rhf       =   B.T0;
B.thf       =   B.T0;
B.bhf       =   B.T0;
% ======================================================================= %
%% ====================== Define time constants ========================= %
T.tmaxini   =   10;          %   Maximum time in Ma
T.itmax     =   100;         %   Maximum number of iterations
T.dtdifac   =   0.8;         %   Diffusion stability criterium
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
Pl.xlab     =   '$$x\ [\ km\ ]$$';
Pl.zlab     =   '$$z\ [\ km\ ]$$';

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
        M.ModDir    = ['data/GaussDiffusion_',...
            '_xmax_',num2str(M.xmax),...
            '_nx_',num2str(N.nz),'_nz_',num2str(N.nz)];
        if ~exist(M.ModDir,'dir')
            mkdir(M.ModDir)
        end
        Pl.filename    = [M.ModDir,'/Evolution.gif'];
        set(figure(1),'position',[1.8,1.8,766.4,780.8]);
        h           =   figure(1);
end
% ======================================================================= %
%% ================ Information for the command window ================== %
fprintf([' Gaussian Diffusion \n  --------------------- ',...
    '\n Diffusion mit: %s',...
    '\n Anfangstemperaturfeld: %s',...
    '\n Aufloesung (nx x nz): %i x %i',...
    '\n  --------------------- ',...
    '\n '],B.DiffMethod,B.Tini,...
    N.nx,N.nz);
fprintf(['Maximum Time : %1.4g',...
    '\n Maximale Anzahl an Iterationen: %5i',...
    '\n  --------------------- \n'],...
    T.tmax,T.itmax)
% ======================================================================= %
T.dtdiff    =   T.dtdifac*min([N.dx,abs(N.dz)])^2/Py.kappa/4;
T.tscale    =   (365.25*3600*24)*1e6;
M.Lscale    =   1e3;
T.dt        =   T.dtdiff;
%% ========================= Time loop ================================= %%
for it = 1:T.itmax
    %% ===================== Calculate time stepping ==================== %
    if it>1
        T.time(it)  =   T.time(it-1) + T.dt;
        
        if T.time(it) > T.tmax
            T.dt        =   T.tmax - T.time(it-1);
            T.time(it)  =   T.time(it-1) + T.dt;
        end
        
        % Analytical Solution ------------------------                    %
        D.Tana      =  B.T0 + B.TAmpl ./...
            (1 + 2 * pi * T.time(it) * Py.kappa / B.Tsigma^2 ) .*...
            exp( -( (M.X-0.5*M.L).^2 + (M.Z-0.5*M.H).^2 ) ./...
            ( 2*B.Tsigma^2 / pi + 4 * T.time(it) * Py.kappa ));
    end
    % =================================================================== %
    %% ========================== Diffusion ============================= %
    if it > 1
        D.T             =   Diffusion(B,D.T,D.Q,D.rho,Py,T.dt,N);
    end
    D.Tmax(it)          =   max(max(D.T));
    D.Tmaxa(it)         =   max(max(D.Tana));
    D.TProfile(:,it)    =   D.T(M.X == M.L/2);
    D.TProfilea(:,it)   =   D.Tana(M.X == M.L/2);
    
    D.epsT      =   (D.Tana - D.T);
    D.RMS(it)   =   sqrt(sum(sum((D.epsT).^2))/(N.nx*N.nz));
    %% ========================== Plot data ============================= %
    Pl.time     =   ...
        [{['@ Iteration: ',sprintf('%i',it)]};...
        {['Time: ',sprintf('%2.2e',T.time(it)/T.tscale),' [Ma]']}];
    
    if (mod(it,5)==0||it==1||T.time(it) >= T.tmax)
        switch Pl.plotfields
            case 'yes'
                figure(1)
                clf
                ax1=subplot(2,2,1);
                plotfield(D.T,M.X./M.Lscale,M.Z./M.Lscale,Pl,'contourf',...
                    '$$T$$','contoury',D.Tana);
                colormap(ax1,flipud(Pl.lajolla))
                ax2=subplot(2,2,2);
                plotfield(D.epsT,M.X./M.Lscale,M.Z./M.Lscale,Pl,'pcolor',...
                    '$$\varepsilon_{T}$$');
                colormap(ax2,Pl.vik)
                caxis([-0.4 0.4])
                ax3=subplot(2,2,3);
                plot(D.TProfile(:,it),M.z./M.Lscale,'k','LineWidth',2)
                hold on
                plot(D.TProfilea(:,it),M.z./M.Lscale,'y--')
                hold off
                axis([B.T0 B.T0+B.TAmpl M.H/M.Lscale 0])
                xlabel('$$T_{x=L/2}\ [^oC]$$','Interpreter','latex')
                ylabel('$$Depth\ [km]$$','Interpreter','latex')
                title([{'$$Temperatur\ Profile$$'};...
                    {'$$at\ Distance\ x=L/2$$'}],'Interpreter','latex')
                set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
                    'TickLabelInterpreter','latex')
                ax4=subplot(2,2,4);
                plot(T.time(1:it)/T.tscale,D.RMS(1:it),'k','LineWidth',2)
                ylabel('$$RMS$$','Interpreter','latex')
                xlabel('$$Time\ [Myrs]$$','Interpreter','latex')
                set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
                    'TickLabelInterpreter','latex')
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
    if (T.time(it) >= T.tmax)
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
% ======================================================================= %
T.tend      = toc(T.tstart);
%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\..\DiffusionProblem')
    rmpath('..\..\SetUp')
else
    rmpath('../../DiffusionProblem')
    rmpath('../../SetUp')
end
% ======================================================================= %
% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


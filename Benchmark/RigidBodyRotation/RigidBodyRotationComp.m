clear
clc
% profile on
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\AdvectionProblem')
    addpath('..\..\SetUp')
    addpath('..\..\StokesProblem')
else
    addpath('../../AdvectionProblem')
    addpath('../../SetUp')
    addpath('../../StokesProblem')
end
% ======================================================================= %
scheme          =   {'upwind','SLF','semi-lag','tracers'};
for i = 1:size(scheme,2)
    %% ===================== Some initial definitions =================== %
    Pl.savefig      =   'yes';
    Pl.plotfields   =   'no';
    % =================================================================== %
    %% ============ Define method to solve the energy equation ========== %
    B.AdvMethod     =   scheme{i};
    B.Aparam        =   'temp';
    % =================================================================== %
    %% ==================== Define viscosity conditions ================= %
    Py.eparam       =   'none';
    B.EtaIni        =   'none';
    % =================================================================== %
    %% ================== Define initial temperature anomaly ============ %
    B.Tini          =   'block';       % Tanomaly
    B.T0            =   1000;           % [ K ]
    B.TAmpl         =   200;            % [ K ]
    switch B.Tini
        case 'gaussianRBR'
            B.Tsigma        =   12;
        otherwise
            B.Tsigma        =   0.1;
    end
    % =================================================================== %
    %% ========================= Define flow field ====================== %
    B.IniFlow       =   'RigidBody';    % FlowField
    B.FlowFac       =   1;
    % ============================================= ====================== %
    %% ==================== Define model geometry constants ============= %
    M.H             =   -1;             % [ km ]
    M.xmax          =   1;              % Aspect ratio
    % =================================================================== %
    %% ====================== Define the numerical grid ================= %
    N.nx            =   101;
    N.nz            =   101;
    % =================================================================== %
    %% ====================== Define time constants ===================== %
    T.tmaxini       =   6.2869e-1;      % [ Ma ]
    T.itmax         =   1e4;
    T.dtfac         =   1.0;            % Courant time factor
    % =================================================================== %
    %% ====================== Tracer advection method =================== %
    N.nmx           =   5;
    N.nmz           =   5;
    % =================================================================== %
    %% ========================= Define fields required ================= %
    [Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
    D                       =   rmfield(D,{'Q','rho','P','Nus','eta'});
    % =================================================================== %
    %% ======================== Setup initial conditions ================ %
    [T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
    % =================================================================== %
    %% =========== Interpolate velocity onto the regular grid =========== %
    [ID]                    =   InterpStaggered(D,ID,N,'velocity');
    % =================================================================== %
    %% ===================== Calculate time stepping ==================== %
    T.dt      =     T.dtfac*min(N.dx,abs(N.dz))/...
        (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));
    %% ========================= Plot parameter ========================= %
    Pl.inc      =   min(N.nz/11,N.nx/11);
    Pl.inc      =   round(Pl.inc);
    Pl.xlab     =   '$$x\ [km]$$';
    Pl.zlab     =   '$$z\ [km]$$';
    % Animation settings ------------------------------------------------ %
    switch Pl.savefig
        case 'yes'
            M.ModDir    =   ['data/',B.AdvMethod,'/'];
            if ~exist(M.ModDir,'dir')
                mkdir(M.ModDir)
            end
            Pl.filename    = [M.ModDir,'2D_Advection_',B.IniFlow,'.gif'];
            if strcmp(B.AdvMethod,'tracers')
                if strcmp(getenv('OS'),'Windows_NT')
                    set(figure(2),'position',[193.0,165.8,1274.4,596.2]);
                else
                    set(figure(2),'position',[170,221,1558,684]);
                end
            else
                set(figure(2),'position',[74.6,198.6,640.8,514.4]);
            end
            h           = figure(2);
    end
    %% ================ Information for the command window ============== %
    fprintf([' Rigid Body Rotation  --------------------- ',...
        '\n Advektion mit: %s',...
        '\n Viskositaet ist: %s',...
        '\n Aufloesung (nx x nz): %i x %i',...
        '\n  --------------------- ',...
        '\n\n '],B.AdvMethod,Py.eparam,...
        N.nx,N.nz);
    fprintf(['Maximum Time : %1.4g',...
        '\n Maximale Anzahl an Iterationen: %5i',...
        '\n  --------------------- \n'],...
        T.tmax/1e6/(365.25*24*60*60),T.itmax)
    % =================================================================== %
    %% ========================= Time loop ============================= %%
    for it = 1:T.itmax
        %         disp([' Time step: ',num2str(it)])
        if it>1
            T.time(it)  =   T.time(it-1) + T.dt;
            if T.time(it) > T.tmax
                T.dt        =   T.tmax - T.time(it-1);
                T.time(it)  =   T.time(it-1) + T.dt;
            end
        end
        %% ========================== Plot data ========================= %
        Pl.time     =   ...
            {};
        switch Pl.plotfields
            case 'yes'
                if (mod(it,50)==0||it==1)
                    switch B.AdvMethod
                        case 'tracers'
                            figure(2),clf
                            tit = {['$$2-D\ numerical\ Advection:',B.AdvMethod,'$$'];...
                                ['$$\Delta\ t_{fac} = ',num2str(T.dtfac),...
                                '; nx = ',num2str(N.nx),', nz = ',...
                                num2str(N.nz),', mpe: ',num2str(N.nmx*N.nmz),...
                                ', Step: ',num2str(it),'$$']};
                            ax1=subplot(1,2,1);
                            plotfield(D.T./D.Tmax,M.X/1e3,M.Z/1e3,Pl,...
                                'pcolor',tit,'quiver',ID.vx,ID.vz)
                            colormap(ax1,flipud(Pl.lajolla))
                            
                            ax2=subplot(1,2,2);
                            plot(Ma.XM./1e3,Ma.ZM./1e3,'.','MarkerSize',1)
                            hold on
                            plot(M.X./1e3,M.Z./1e3,'kx','MarkerSize',2)
                            xlabel('$$x\ [km]$$','Interpreter','latex')
                            ylabel('$$z\ [km]$$','Interpreter','latex')
                            title('$$Tracerdistribution$$','Interpreter','latex')
                            axis equal; axis tight
                            set(gca,'FontWeight','Bold','TickLabelInterpreter','latex')
                        otherwise
                            figure(2),clf
                            tit = {['$$2-D\ numerical\ Advection:',B.AdvMethod,'$$'];...
                                ['$$\Delta\ t_{fac} = ',num2str(T.dtfac),...
                                '; nx = ',num2str(N.nx),', nz = ',...
                                num2str(N.nz),', mpe: ',num2str(N.nmx*N.nmz),...
                                ', Step: ',num2str(it),'$$']};
                            ax1=subplot(1,1,1);
                            plotfield(D.T./D.Tmax,M.X/1e3,M.Z/1e3,Pl,...
                                'pcolor',tit,'quiver',ID.vx,ID.vz)
                            caxis([B.T0/(B.T0+B.TAmpl) 1])
                            colormap(ax1,flipud(Pl.lajolla))
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
                end
            case 'no'
                disp(['Iteration: ',sprintf('%i',it)])
        end
        %% ========================== Advection ========================= %
        [D,Ma,ID]       =   Advection(it,N,B,D,ID,Py,T.dt,M,Ma);
        % =============================================================== %
        if T.time(it) >= T.tmax
            switch B.AdvMethod
                case 'tracers'
                    figure(2),clf
                    tit = {['$$2-D\ numerical\ Advection:',B.AdvMethod,'$$'];...
                        ['$$\Delta\ t_{fac} = ',num2str(T.dtfac),...
                        '; nx = ',num2str(N.nx),', nz = ',...
                        num2str(N.nz),', mpe: ',num2str(N.nmx*N.nmz),...
                        ', Step: ',num2str(it),'$$']};
                    ax1=subplot(1,2,1);
                    plotfield(D.T./D.Tmax,M.X/1e3,M.Z/1e3,Pl,...
                        'pcolor',tit,'quiver',ID.vx,ID.vz)
                    colormap(ax1,flipud(Pl.lajolla))
                    
                    ax2=subplot(1,2,2);
                    plot(Ma.XM./1e3,Ma.ZM./1e3,'.','MarkerSize',1)
                    hold on
                    plot(M.X./1e3,M.Z./1e3,'kx','MarkerSize',2)
                    xlabel('$$x\ [km]$$','Interpreter','latex')
                    ylabel('$$z\ [km]$$','Interpreter','latex')
                    title('$$Tracerdistribution$$','Interpreter','latex')
                    axis equal; axis tight
                    set(gca,'FontWeight','Bold','TickLabelInterpreter','latex')
                otherwise
                    figure(2),clf
                    tit = {['$$2-D\ numerical\ Advection:',B.AdvMethod,'$$'];...
                        ['$$\Delta\ t_{fac} = ',num2str(T.dtfac),...
                        '; nx = ',num2str(N.nx),', nz = ',...
                        num2str(N.nz),', mpe: ',num2str(N.nmx*N.nmz),...
                        ', Step: ',num2str(it),'$$']};
                    ax1=subplot(1,1,1);
                    plotfield(D.T./D.Tmax,M.X/1e3,M.Z/1e3,Pl,...
                        'pcolor',tit,'quiver',ID.vx,ID.vz)
                    caxis([B.T0/(B.T0+B.TAmpl) 1])
                    colormap(ax1,flipud(Pl.lajolla))
            end
            break
        end
    end
    figure(3)
    tit = {['$$2-D\ numerical\ Advection:',B.AdvMethod,'$$'];...
        ['$$\Delta\ t_{fac} = ',num2str(T.dtfac),...
        '; nx = ',num2str(N.nx),', nz = ',...
        num2str(N.nz),', mpe: ',num2str(N.nmx*N.nmz),...
        ', Step: ',num2str(it),'$$']};
    ax1=subplot(2,2,i);
    plotfield(D.T./D.Tmax,M.X/1e3,M.Z/1e3,Pl,...
        'pcolor',tit,'quiver',ID.vx,ID.vz)
    caxis([B.T0/(B.T0+B.TAmpl) 1])
    colormap(ax1,flipud(Pl.lajolla))
end
%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\..\AdvectionProblem')
    rmpath('..\..\SetUp')
    rmpath('..\..\StokesProblem')
else
    rmpath('../../AdvectionProblem')
    rmpath('../../SetUp')
    rmpath('../../StokesProblem')
end
% ======================================================================= %
% profile viewer
% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %

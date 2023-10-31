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

Schema          =   {'explicit';'implicit';'CNV';'ADI'};
nrnxnz          =   2;

eps             =   zeros(size(Schema,1),nrnxnz);
nxnz            =   zeros(size(Schema,1),nrnxnz);
Tmax            =   zeros(size(Schema,1),nrnxnz);
Tmean           =   zeros(size(Schema,1),nrnxnz);

for k = 1:size(Schema,1)
    disp('')
    disp(Schema{k})
    for l = 1:nrnxnz
        
        %% ===================== Some initial definitions ======================= %
        Pl.savefig      =   'no';
        Pl.plotfields   =   'yes';
        Py.scale        =   'no';
        % ======================================================================= %
        %% ============ Define method to solve the energy equation ============== %
        B.DiffMethod    =   Schema{k};
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
        N.nz        =   l*20+1;             %   Vertical grid solution
        N.nx        =   l*20+1;             %   Horizontal grid solution
        disp(['nx = ',num2str(N.nx)])
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
        B.lhf       =   B.T0;   B.rhf       =   B.T0;
        B.thf       =   B.T0;   B.bhf       =   B.T0;
        % ======================================================================= %
        %% ====================== Define time constants ========================= %
        T.tmaxini   =   10;           %   Maximum time in Ma
        T.itmax     =   5000;        %   Maximum number of iterations
        T.dtdifac   =   0.95;         %   Diffusion stability criterium
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
        Pl.xlab     =   '$$x\ [km]$$';
        Pl.zlab     =   '$$z\ [km]$$';
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
                M.ModDir    = ['data/ResTest/GaussDiffusion_',...
                    B.DiffMethod,'_xmax_',num2str(M.xmax),...
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
            D.Tmean(it)         =   mean(mean(D.T));
            D.Tmaxa(it)         =   max(max(D.Tana));
            D.TProfile(:,it)    =   D.T(M.X == M.L/2);
            D.TProfilea(:,it)   =   D.Tana(M.X == M.L/2);
            
            D.epsT      =   (D.Tana - D.T);
            D.RMS(it)   =   sqrt(sum(sum((D.epsT).^2))/(N.nx*N.nz));
            %% ========================== Plot data ============================= %
            Pl.time     =   ...
                [{['@ Iteration: ',sprintf('%i',it)]};...
                {['Time: ',sprintf('%2.2e',T.time(it)/T.tscale),' [Ma]']}];
            
            if (mod(it,2)==0||it==1||T.time(it)==T.tmax)
                switch Pl.plotfields
                    case 'yes'
                        figure(1)
                        clf
                        ax1=subplot(2,2,1);
                        plotfield(D.T,M.X./M.Lscale,M.Z./M.Lscale,Pl,'contourf',...
                            '$$T$$','contoury',D.Tana);
                        colormap(ax1,Pl.lajolla)
                        ax2=subplot(2,2,2);
                        plotfield(D.epsT,M.X./M.Lscale,M.Z./M.Lscale,Pl,'pcolor',...
                            '$$\varepsilon_{T}$$');
                        colormap(ax2,Pl.vik)
                        caxis([-0.5 0.5])
                        ax3=subplot(2,2,3);
                        plot(D.TProfile(:,it),M.z./M.Lscale,'k','LineWidth',2)
                        hold on
                        plot(D.TProfilea(:,it),M.z./M.Lscale,'y--')
                        hold off
                        axis([B.T0 B.T0+B.TAmpl M.H/M.Lscale 0])
                        xlabel('$$T_{x=L/2}\ [^oC]$$','Interpreter','latex')
                        ylabel('$$Depth\ [km]$$','Interpreter','latex')
                        title([{'$$Temperature\ Profile$$'};...
                            {'$$@\ Distance\ x=L/2$$'}],'Interpreter','latex')
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
            % =================================================================== %
            if T.time(it) >= T.tmax
                T.itmax = it;
                break
            end
        end
        %% ======================== Save final figure =========================== %
        switch Pl.savefig
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
                caxis([-0.5 0.5])
                ax3=subplot(2,2,3);
                plot(D.TProfile(:,it),M.z./M.Lscale,'k','LineWidth',2)
                hold on
                plot(D.TProfilea(:,it),M.z./M.Lscale,'y--')
                hold off
                axis([B.T0 B.T0+B.TAmpl M.H/M.Lscale 0])
                xlabel('$$T_{x=L/2}\ [^oC]$$','Interpreter','latex')
                ylabel('$$Depth\ [km]$$','Interpreter','latex')
                title([{'$$Temperature\ Profile$$'};...
                    {'$$@\ Distance\ x=L/2$$'}],'Interpreter','latex')
                set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
                    'TickLabelInterpreter','latex')
                ax4=subplot(2,2,4);
                plot(T.time(1:it)/T.tscale,D.RMS(1:it),'k','LineWidth',2)
                ylabel('$$RMS$$','Interpreter','latex')
                xlabel('$$Time\ [Myrs]$$','Interpreter','latex')
                set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
                    'TickLabelInterpreter','latex')
                saveas(figure(1),[M.ModDir,'/Field_SS'],'png')
        end
        % ======================================================================= %
        disp(' Time loop finished ... ')
        disp('-> Use new grid size...')
        
        eps(k,l)    = D.RMS(T.itmax);
        nxnz(k,l)   = N.nx*N.nz;
        Tmax(k,l)   = D.Tmax(T.itmax);
        Tmean(k,l)  = D.Tmean(T.itmax);
    end
end
T.tend      = toc(T.tstart);

legendinfo  = cell(1,size(Schema,1)+1);
linstyle    = {'-','--',':','-.','-'};
figure(2)
clf
for k = 1:size(Schema,1)
    subplot(1,3,1)
    p(k) = loglog(1./nxnz(k,:),eps(k,:),'LineStyle',linstyle{k},...
        'Marker','*','LineWidth',2);
    legendinfo{k} = ['$$',Schema{k},'$$'];
    hold on
    if k == size(Schema,1)
        legendinfo{k+1} = '';
        xlabel('$$1/nx/nz$$','Interpreter','latex')
        ylabel('$$RMS_{\Delta T}$$','Interpreter','latex')
        axis square
        set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
            'TickLabelInterpreter','latex')
        legend(p,legendinfo,'Location','NorthWest','Interpreter','latex')
        %         axis([1e-5 1e-2 1e-2 2])
    end
    subplot(1,3,2)
    loglog(1./nxnz(k,:),Tmax(k,:),'LineStyle',linstyle{k},...
        'Marker','*','LineWidth',2)
    hold on
    if k == size(Schema,1)
        legendinfo{k+1} = '$$Sol_{ana}$$';
        plot(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)),...
            max(max(D.Tana)).*...
            ones(1,length(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)))),'k--')
        xlabel('$$1/nx/nz$$','Interpreter','latex')
        ylabel('$$T_{max}$$','Interpreter','latex')
        axis square
        set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
            'xscale','log','TickLabelInterpreter','latex')
        legend(legendinfo,'Location','NorthWest','Interpreter','latex')
        %         axis([1e-5 1e-2 785 795])
    end
    subplot(1,3,3)
    loglog(1./nxnz(k,:),Tmean(k,:),'LineStyle',linstyle{k},...
        'Marker','*','LineWidth',2)
    hold on
    if k == size(Schema,1)
        plot(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)),...
            mean(mean(D.Tana)).*...
            ones(1,length(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)))),'k--')
        xlabel('$$1/nx/nz$$','Interpreter','latex')
        ylabel('$$\langle T \rangle$$','Interpreter','latex')
        axis square
        set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2,...
            'TickLabelInterpreter','latex')
        legend(legendinfo,'Location','NorthEast','Interpreter','latex')
        %         axis([1e-5 1e-2 775 776])
    end
end

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


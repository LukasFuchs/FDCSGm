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
M.H             =   -1e3;                 %   Modeltiefe [ in m ]
% ======================================================================= %
%% Parameters to calculate the growth rate ============================== %
lambda          =   [0.5 0.75 1.0 1.25 1.5 1.75 ...
    2.0 2.25 2.5 2.75 ...
    3.0 3.25 3.5 3.75 ...
    4.0].*1e3;              % [ m ]
Py.eta1         =   1e18;   % Viscosity lower layer [ Pa s ]
etar            =   [1e-6 1 10 100 500];    % eta1/eta0
b1              =   [0.5 1 5 50 250];
b2              =   [0.2 0.15 0.1 0.05 0];

PP.K            =   zeros(length(lambda),length(etar));
PP.phi          =   zeros(length(lambda),length(etar));

%% Analytical solution ================================================== %
lambda_ana      =   linspace(0.3,5,51).*1e3;
PP.K_ana        =   zeros(length(lambda_ana),length(etar));
% ======================================================================= %
for i = 1:length(etar)
    Py.eta0             =   Py.eta1*etar(i);   %   Viskositaet comp. 0 [ Pa*s ]
    phi1     =   (2.*pi.*(abs(M.H)/2))./lambda_ana;
    phi2     =   (2.*pi.*(abs(M.H)/2))./lambda_ana;
    c11         =   (Py.eta0.*2.*phi1.^2)./...
        (Py.eta1.*(cosh(2.*phi1) - 1 - 2.*phi1.^2)) - ...
        (2.*phi2.^2)./...
        (cosh(2.*phi2) - 1 - 2.*phi2.^2);
    d12     =   (Py.eta0.*(sinh(2.*phi1) - 2.*phi1))./...
        (Py.eta1.*(cosh(2.*phi1) - 1 - 2.*phi1.^2)) + ...
        (sinh(2.*phi2) - 2.*phi2)./...
        (cosh(2.*phi2) - 1 - 2.*phi2.^2);
    i21     =   (Py.eta0.*phi2.*(sinh(2.*phi1) + 2.*phi1))./...
        (Py.eta1.*(cosh(2.*phi1) - 1 - 2.*phi1.^2)) + ...
        (phi2.*(sinh(2.*phi2) + 2.*phi2))./...
        (cosh(2.*phi2) - 1 - 2.*phi2.^2);
    j22     =   (Py.eta0.*2.*phi1.^2.*phi2)./...
        (Py.eta1.*(cosh(2.*phi1) - 1 - 2.*phi1.^2)) - ...
        (2.*phi2.^3)./...
        (cosh(2.*phi2) - 1 - 2.*phi2.^2);
    
    PP.K_ana(:,i)   =   -d12./(c11.*j22-d12.*i21);
    PP.phi_ana      =   phi1;
    for j = 1:length(lambda)
        M.H             =   M.H/1e3;            % in [km]
        %% ==================== Define viscosity conditions ===================== %
        Py.eparam       =   'variable';
        B.EtaIni        =   'RTI';
        B.lambda        =   lambda(j)/1e3;      % Wellenlaenge der Perturbation [ km ]
        % ======================================================================= %
        %% ========================= Define flow field ========================== %
        B.FlowFac       =   10;
        % ======================================================================= %
        %% ==================== Define model geometry constants ================= %
        M.xmax          =   -2*B.lambda/M.H;    %   Seitenverhaeltniss
        B.deltaA        =   -M.H*1e3/1500;        % Amplitude [ m ]
        % ======================================================================= %
        %% ====================== Define the numerical grid ===================== %
        N.nz            =   101;                 %   Vertikale Gitteraufloesung
        N.nx            =   (N.nz-1)*M.xmax+1;   %   Horizontale Gitteraufloesung
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
        T.itmax     =   1;          %   Maximal erlaubte Anzahl der Iterationen
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
                    set(figure(1),'position',[1.0,1.0,1536.0,788.8]);
                    h           =   figure(1);
                else
                    set(figure(1),'position',[-1919,1,960,988]);
                    h           =   figure(1);
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
        %     for it = 1:T.itmax
        if(strcmp(B.AdvMethod,'none')==0)
            [D,A]       =   solveSECE(D,Py,N,B);
        end
        % =================================================================== %
        %     %% ===================== Calculate time stepping ==================== %
        %     T.dt        =   T.dtfac*min(N.dx,abs(N.dz))/...
        %         (sqrt(max(max(D.vx))^2 + max(max(D.vz))^2));
        %     if it>1
        %         T.time(it)  =   T.time(it-1) + T.dt;
        %     end
        %     % =================================================================== %
        %% =========== Interpolate velocity onto the regular grid =========== %
        [ID]        =   InterpStaggered(D,ID,N,'velocity');
        %         D.meanV(it) =   rms(ID.vx(:) + ID.vz(:));
        % =================================================================== %
        %% ========================== Plot data ============================= %
        Pl.time     =   ...
            [];
        switch Pl.plotfields
            case 'yes'
                %             if (mod(it,5)==0||it==1)
                figure(1)
                clf
                ax1=subplot(1,1,1);
                %                 Pl.cbtitle  =   [{'$$\rho$$'},{'$$[kg/m^3]$$'}];
                %                 plotfield(D.rho,M.X/1e3,M.Z/1e3,Pl,'pcolor',...
                %                     [],'quiver',ID.vx,ID.vz)
                %                 colormap(ax1,flipud(Pl.oslo))
                ax2=subplot(1,1,1);
                Pl.cbtitle  =   [{'$$\eta$$'},{'$$[ Pa s ]$$'}];
                plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
                    [],'quiver',ID.vx,ID.vz)
                colormap(ax2,flipud(Pl.lapaz))
                %                 ax3=subplot(3,1,3);
                %                 Pl.cbtitle  =   [{'v'},{'[ cm/a ]'}];
                %                 plotfield(ID.v.*100.*(60*60*24*365.25),...
                %                     M.X/1e3,M.Z/1e3,Pl,'pcolor',[])
                %                 colormap(ax3,Pl.imola)
                
                switch Pl.savefig
                    case 'yes'
                        saveas(figure(1),...
                            [M.ModDir,'/Field',num2str(it)],'png')
                        %                         % Capture the plot as an image
                        %                         frame       = getframe(h);
                        %                         im          = frame2im(frame);
                        %                         [imind,cm]  = rgb2ind(im,256);
                        %
                        %                         % Write to the GIF File
                        %                         if it == 1
                        %                             imwrite(imind,cm,Pl.filename,'gif', 'Loopcount',inf);
                        %                         else
                        %                             imwrite(imind,cm,Pl.filename,'gif','WriteMode','append');
                        %                         end
                end
                %             end
            case 'no'
                disp(['Iteration: ',sprintf('%i',it)])
        end
        % =================================================================== %
        %     %% ========================== Advection ============================= %
        %     [D,Ma,ID]       =   Advection(it,N,B,D,ID,Py,T.dt,M,Ma);
        %     % =================================================================== %
        %     end
        
        % Position
        PP.xwave        =   M.xmax*1e3/2;       % [ m ]
        PP.zwave        =   M.H/2 + B.deltaA;   % [ m ]
        
        PP.xn           =   double(int16(PP.xwave./N.dx));
        PP.zn           =   double(int16(PP.zwave./N.dz))+1;
        
        PP.dx           =   (PP.xwave+N.dx/2)./N.dx-PP.xn;
        PP.dz           =   PP.zwave./N.dz-PP.zn+1;
        
        PP.wvy          =   (1.0-PP.dx)*(1.0-PP.dz)*D.vz(PP.zn,PP.xn) + ...
            PP.dx*(1.0-PP.dz)*D.vz(PP.zn,PP.xn+1) + ...
            (1.0-PP.dx)*PP.dz*D.vz(PP.zn+1,PP.xn) + ...
            PP.dx*PP.dz*D.vz(PP.zn+1,PP.xn+1);
        %     hold on
        %     plot(PP.xwave/1e3,PP.zwave/1e3,'o','MarkerSize',8)
        %     plot(M.x1(PP.xn)/1e3,M.z(PP.zn)/1e3,'d','MarkerSize',8)
        %     plot(M.X1(PP.zn,PP.xn)/1e3,M.Z(PP.zn,PP.xn)/1e3,'r.',...
        %         M.X1(PP.zn,PP.xn+1)/1e3,M.Z(PP.zn,PP.xn+1)/1e3,'r.',...
        %         M.X1(PP.zn+1,PP.xn)/1e3,M.Z(PP.zn+1,PP.xn)/1e3,'r.',...
        %         M.X1(PP.zn+1,PP.xn+1)/1e3,M.Z(PP.zn+1,PP.xn+1)/1e3,'r.')
        PP.Q        =   (Py.rho0-Py.rho1)*(abs(M.H))/2*Py.g/2/Py.eta1;
        PP.K(j,i)   =   PP.wvy/B.deltaA/PP.Q;
        PP.phi(j,i) =   2*pi*(abs(M.H)/2)/(B.lambda*1e3);
        
    end
    figure(2)
    plot(PP.phi(:,i),b1(i).*PP.K(:,i) + b2(i),'ko','MarkerSize',8)
    hold on
    plot(PP.phi_ana,b1(i).*PP.K_ana(:,i) + b2(i),'-')
    xlabel('$$\phi_1 = 2\pi h_1/\lambda$$','Interpreter','latex')
    ylabel('$$b_1K+b_2$$','Interpreter','latex')
    title('','Interpreter','latex')
    axis([0.5 4 0.05 0.4])
    set(gca,'FontWeight','Bold','FontSize',12,...
    'TickLabelInterpreter','latex')
end
% %% ======================== Save final figure =========================== %
% figure(3)
% clf
% ax1=subplot(3,1,1);
% Pl.cbtitle  =   [{'$$\rho$$'},{'$$[kg/m^3]$$'}];
% plotfield(D.rho,M.X/1e3,M.Z/1e3,Pl,'contourf',...
%     [],'quiver',ID.vx,ID.vz)
% colormap(ax1,flipud(Pl.oslo))
% ax2=subplot(3,1,2);
% Pl.cbtitle  =   [{'$$\eta$$'},{'$$[ Pa s ]$$'}];
% plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
%     [])
% colormap(ax2,flipud(Pl.lapaz))
% ax3=subplot(3,1,3);
% Pl.cbtitle  =   [{'v'},{'[ cm/a ]'}];
% plotfield(ID.v.*100.*(60*60*24*365.25),...
%     M.X/1e3,M.Z/1e3,Pl,'pcolor',[])
% colormap(ax3,Pl.imola)
% switch Pl.savefig
%     case 'yes'
%         saveas(figure(3),['data/RTI_adv_',B.AdvMethod,...
%             '_eta_',Py.eparam,'_xmax_',num2str(M.xmax),...
%             '_nz_',num2str(N.nz),'_SS'],'png')
% end
% % ======================================================================= %
% %% ======================== Plot time serieses ========================== %
% figure(4)
% plot(T.time/1e6/(365.25*24*60*60),...
%     D.meanV*100*365.25*24*60*60,...
%     'LineWidth',2)
% set(gca,'FontWeight','Bold',...
%     'LineWidth',2,'FontSize',15,'TickLabelInterpreter','latex')
% xlabel('$$t\ [Ma]$$','Interpreter','latex')
% ylabel('$$V_{RMS}\ [cm/a]$$','Interpreter','latex')
% title('$$Root\ Mean\ Square\ Velocity$$','Interpreter','latex')
%
% switch Pl.savefig
%     case 'yes'
%         saveas(figure(4),[M.ModDir,'/TimeSeries'],'png')
% end
% % ======================================================================= %
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


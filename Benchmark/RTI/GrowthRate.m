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
Pl.plotfields   =   'no';
% ======================================================================= %
%% ============ Define method to solve the energy equation ============== %
B.AdvMethod     =   'tracers';
B.Aparam        =   'comp';
M.H             =   -1e3;                 %   Modeltiefe [ in m ]
% ======================================================================= %
%% Parameters to calculate the growth rate ============================== %
% Wavelength for the compositional perturabation ------------------------ %
lambda          =   [0.5 0.75 1.0 1.25 1.5 1.75 ...
    2.0 2.25 2.5 2.75 ...
    3.0 3.25 3.5 3.75 ...
    4.0].*1e3;              % [ m ]
% Viscosity of the lower layer ------------------------------------------ %
Py.eta1         =   1e18;                       % [ Pa s ]
etar            =   [1e-6 1 10 100 500];        % eta1/eta0
% Plotting factors following Gerya (2009) ------------------------------- %
b1              =   [0.5 1 5 50 250];
b2              =   [0.2 0.15 0.1 0.05 0];
% Divisional factor of the amplitude following Gerya (2009) ------------- %
delfac          =   [15 150];
% Plotting modifications ------------------------------------------------ %
mark            =   {'o','.'};
marksz          =   {8,10};
% Array definition ------------------------------------------------------ %
PP.K            =   zeros(length(lambda),length(etar));
PP.phi          =   zeros(length(lambda),length(etar));
%% Analytical solution ================================================== %
lambda_ana      =   linspace(0.3,5,51).*1e3;
PP.K_ana        =   zeros(length(lambda_ana),length(etar));
% ======================================================================= %
for k = 1:length(delfac)
    for i = 1:length(etar)
        Py.eta0             =   Py.eta1*etar(i);
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
            %% ==================== Define viscosity conditions ========= %
            Py.eparam       =   'variable';
            B.EtaIni        =   'RTI';
            B.lambda        =   lambda(j)/1e3;  % [ km ]
            % =========================================================== %
            %% ==================== Define flow field =================== %
            B.FlowFac       =   10;
            % =========================================================== %
            %% ============= Define model geometry constants ============ %
            M.xmax          =   -2*B.lambda/M.H;
            B.deltaA        =   -(M.H/2)*1e3/delfac(k);
            % =========================================================== %
            %% =============== Define the numerical grid ================ %
            N.nz            =   51;
            N.nx            =   (N.nz-1)*M.xmax+1;
            % =========================================================== %
            %% ================ Tracer advection method ================= %
            % So far only required for initialization!
            N.nmx       =   5;
            N.nmz       =   5;
            % =========================================================== %
            %% ================ Define physical constants =============== %
            Py.g        =   10.0;           %   Schwerebeschleunigung [m/s^2]
            Py.rho0     =   3300;           %   Hintergunddichte [kg/m^3]
            Py.k        =   3;              %   Thermische Leitfaehigkeit [ W/m/K ]
            Py.cp       =   1000;           %   Heat capacity [ J/kg/K ]
            Py.alpha    =   5e-5;           %   Thermischer Expnasionskoef. [ K^-1 ]
            
            Py.kappa    =   Py.k/Py.rho0/Py.cp; % 	Thermische Diffusivitaet [ m^2/s ]
            
            Py.Q0       =   0;              %   Waermeproduktionsrate pro Volumen [W/m^3]
            
            Py.rho1     =   3000;           %   Density comp. 1 [ kg/m^3 ]
            % =========================================================== %
            %% =============== Define boundary conditions =============== %
            %   0 - no slip; 1 - free slip;
            B.tbc       =   0;          %   Top
            B.bbc       =   0;          %   Bottom
            B.lbc       =   1;          %   Left
            B.rbc       =   1;          %   Right
            % =========================================================== %
            %% ====================== Define time constants ============= %
            %             T.tmaxini   =   4500;       %   Maximum time in Ma
            %             T.itmax     =   1;          %   Maximal erlaubte Anzahl der Iterationen
            %             T.dtfac     =   1.0;        %   Advektionscourantkriterium
            T           =   [];
            % =========================================================== %
            %% ==================== Define fields required ============== %
            [Py,D,ID,M,N,T,A,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
            % =========================================================== %
            %% =================== Setup initial conditions ============= %
            [T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
            % =========================================================== %
            %% ========================= Plot parameter ================= %
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
            % =========================================================== %
            %% =========== Information for the command window =========== %
            fprintf([' Rayleigh-Taylor Instabilitaet  -------------- ',...
                '\n Advektion mit: %s',...
                '\n Viskositaet ist: %s',...
                '\n Aufloesung (nx x nz): %i x %i',...
                '\n Refrenzviskositaet [Pa s]: %2.2e',...
                '\n  --------------------- ',...
                '\n\n '],B.AdvMethod,Py.eparam,...
                N.nx,N.nz,Py.eta0);
            % =========================================================== %
            %% ========================= Solving ======================= %%
            if(strcmp(B.AdvMethod,'none')==0)
                [D,A]       =   solveSECE(D,Py,N,B);
            end
            % =========================================================== %
            %% ======== Interpolate velocity onto the regular grid ====== %
            [ID]        =   InterpStaggered(D,ID,N,'velocity');
            % =========================================================== %
            %% ========================== Plot data ===================== %
            Pl.time     =   ...
                [];
            switch Pl.plotfields
                case 'yes'
                    figure(1)
                    clf
                    ax1=subplot(1,1,1);
                    Pl.cbtitle  =   [{'$$\eta$$'},{'$$[ Pa s ]$$'}];
                    plotfield(log10(D.eta),M.X/1e3,M.Z/1e3,Pl,'pcolor',...
                        [],'quiver',ID.vx,ID.vz)
                    colormap(ax1,flipud(Pl.lapaz))
            end
            % =========================================================== %
            % Position
            PP.xwave        =   M.xmax*1e3/2;       % [ m ]
            PP.zwave        =   M.H/2 + B.deltaA;   % [ m ]
            
            PP.xn           =   double(int16(PP.xwave./N.dx));
            PP.zn           =   double(int16(PP.zwave./N.dz));
            
            PP.dx           =   (PP.xwave+N.dx/2)./N.dx-PP.xn;
            PP.dz           =   PP.zwave./N.dz-PP.zn+1;
            
            PP.wvy          =   (1.0-PP.dx)*(1.0-PP.dz)*D.vz(PP.zn,PP.xn) + ...
                PP.dx*(1.0-PP.dz)*D.vz(PP.zn,PP.xn+1) + ...
                (1.0-PP.dx)*PP.dz*D.vz(PP.zn+1,PP.xn) + ...
                PP.dx*PP.dz*D.vz(PP.zn+1,PP.xn+1);
            %         hold on
            %         axis equal; axis tight
            %         plot(PP.xwave/1e3,PP.zwave/1e3,'o','MarkerSize',8)
            %         plot(M.x1(PP.xn)/1e3,M.z(PP.zn)/1e3,'d','MarkerSize',8)
            %         plot(M.X1(PP.zn,PP.xn)/1e3,M.Z(PP.zn,PP.xn)/1e3,'r.',...
            %             M.X1(PP.zn,PP.xn+1)/1e3,M.Z(PP.zn,PP.xn+1)/1e3,'r.',...
            %             M.X1(PP.zn+1,PP.xn)/1e3,M.Z(PP.zn+1,PP.xn)/1e3,'r.',...
            %             M.X1(PP.zn+1,PP.xn+1)/1e3,M.Z(PP.zn+1,PP.xn+1)/1e3,'r.')
            %         ind     =   Ma.XM >= (PP.xwave-100) & ...
            %                     Ma.XM <= (PP.xwave+100) & ...
            %                     Ma.ZM >= (M.H/2+0.1*M.H) & ...
            %                     Ma.ZM <= (M.H/2-0.1*M.H);
            %         scatter(Ma.XM(ind)./1e3,Ma.ZM(ind)./1e3,8,Ma.C(ind),'filled')
            %         axis([(PP.xwave-100)/1e3 (PP.xwave+100)/1e3...
            %             (M.H/2+0.1*M.H)/1e3 (M.H/2-0.1*M.H)/1e3])
            PP.Q            =   (Py.rho0-Py.rho1)*(abs(M.H)/2)*Py.g/2/Py.eta1;
            PP.K(k,j,i)     =   PP.wvy/B.deltaA/PP.Q;
            PP.phi(k,j,i)   =   2*pi*(abs(M.H)/2)/(B.lambda*1e3);
        end
        figure(2)
        plot(PP.phi(k,:,i),b1(i).*PP.K(k,:,i) + b2(i),'LineStyle','none',...
            'Marker',mark{k},'MarkerSize',marksz{k},'MarkerEdgeColor','k')
        hold on
        if k == 1
            plot(PP.phi_ana,b1(i).*PP.K_ana(:,i) + b2(i),'k-')
            xlabel('$$\phi_1 = 2\pi h_1/\lambda$$','Interpreter','latex')
            ylabel('$$b_1K+b_2$$','Interpreter','latex')
            title('','Interpreter','latex')
            axis([0.5 4 0.05 0.4])
            set(gca,'FontWeight','Bold','FontSize',12,...
                'TickLabelInterpreter','latex')
        end
    end
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


clear
clc
%% ======================== Add required paths ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    addpath('..\..\StokesProblem')
    addpath('..\..\SetUp')
else
    addpath('../../StokesProblem')
    addpath('../../SetUp')
    addpath('../../ScaleParam')
end
% ======================================================================= %
%% ===================== Some initial definitions ======================= %
Pl.savefig      =   'no';
Pl.plotfields   =   'yes';
% ======================================================================= %
%% ==================== Define viscosity conditions ===================== %
Py.eparam       =   'variable';
B.EtaIni        =   'exp';
Py.eta0         =   1e20;           % Viscosity at the bottom
Py.eta1         =   1e26;           % Viscosity at the top
% ======================================================================= %
%% ========================= Define flow field ========================== %
B.IniFlow       =   'Channel';
Py.dPdx         =   -200;            % Horizontal pressure gradient; -10.0;
Py.v0           =   1.58e-9;        % Velocity at the top [ m/s ]
B.FlowFac       =   10;
% ----------------------------------------------------------------------- %
%% ==================== Define model geometry constants ================= %
M.H         =   -400;               %   Depth [ in km ]
M.xmax      =   2;                  %   Aspect ratio
% ======================================================================= %
%% ====================== Define the numerical grid ===================== %
N.nz        =   101;             %   Vertical grid solution
N.nx        =   201;             %   Horizontal grid solution
% ======================================================================= %
%% ====================== Define physical constants ===================== %
Py.g        =   10;                 %   Gravitational acceleration [m/s^2]
Py.rho0     =   4000;               %   Reference density [kg/m^3]
% ======================================================================= %
%% ===================== Define boundary conditions ===================== %
% Velocity boundary conditions ------------------------------------------ %
%   0 - no slip; 1 - free slip;
B.tbc       =   0;              %   Top
B.bbc       =   0;              %   Bottom
B.lbc       =   0;              %   Left
B.rbc       =   0;              %   Right
% ======================================================================= %
%% ====================== Define time constants ========================= %
T                       =   [];
% % ======================================================================= %
%% ========================= Define fields required ===================== %
[Py,D,ID,M,N,T,~,Pl]    =   SetUpFields(Py,B,N,M,T,Pl);
% ======================================================================= %
%% ======================== Setup initial conditions ==================== %
[T,D,B,M,Ma,Py]         =   SetUpInitialConditions(T,D,Py,M,N,B);
% ======================================================================= %
%% ========================= Plot parameter ============================= %
Pl.inc      =   min(N.nz/10,N.nx/5);
Pl.inc      =   round(Pl.inc);
Pl.xlab     =   'x [ km ]';
Pl.zlab     =   'z [ km ]';
Pl.tstpinc  =   1;
Pl.cbtitle  =   '$$log_{10} ( \eta )$$';
switch Pl.plotfields
    case 'yes'
        if strcmp(getenv('OS'),'Windows_NT')
            set(figure(1),'position',[1.8,1.8,766.4,780.8]);
            h           =   figure(1);
            clf
        else
            set(figure(1),'position',[-1919,1,960,988]);
            h           =   figure(1);
            clf
        end
end
% Animation settings ---------------------------------------------------- %
switch Pl.savefig
    case 'yes'
        M.ModDir    = 'data/';
        if ~exist(M.ModDir,'dir')
            mkdir(M.ModDir)
        end        
end
% ======================================================================= %
%% ================ Information for the command window ================== %
fprintf(['... \n  --------------------- ',...
    '\n Viskositaet ist: %s',...
    '\n Aufloesung (nx x nz): %i x %i',...
    '\n Refrenzviskositaet [Pa s]: %2.2e',...
    '\n  --------------------- ',...
    '\n '],Py.eparam,...
    N.nx,N.nz,Py.eta0);
% ======================================================================= %
%% =========================== Solving ================================= %%
[D,A]       =   solveSECE(D,Py,N,B);
D.vx_mean   =   mean(D.vx(1:N.nz1,:),2);
% ======================================================================= %
%% =============== Interpolate velocity onto the regular grid =========== %
[ID]        =   InterpStaggered(D,ID,N,'velocity');
% ======================================================================= %
%% ========================= Analytical solution ======================== %
m           =   Py.eta1/Py.eta0;    % eta_top/eta_bottom
if m == 1
    D.vx_ana    =   -Py.dPdx/2/Py.eta0*(M.H.*M.z1' - M.z1'.^2) + ...
        Py.v0.*M.z1'./M.H;
else
    D.vx_ana    =   -Py.dPdx*M.H/(Py.eta0)/log(m).*...
        (m.^(-M.z1'./M.H)/(m-1).*(M.z1'.*(m-1) + M.H) - M.H/(m-1)) - ...
        m.^(-M.z1'./M.H).*m.*Py.v0./(m-1) + ...
        Py.v0*m/(m-1);
end
D.vx_ana        =   flipud(D.vx_ana);
D.Deltav        =   sqrt((D.vx_mean-D.vx_ana).^2./max(D.vx_ana).^2);
% ======================================================================= %
%% ========================== Plot data ================================= %
Pl.time     =   [];
ax1=subplot(2,1,1);
plotfield(log10(D.eta),M.X./1e3,M.Z./1e3,Pl,'contourf',...
    'v','quiver',ID.vx,ID.vz)
colormap(ax1,flipud(Pl.lapaz))
ax2=subplot(2,2,3);
py  =   [M.z1./1e3,fliplr(M.z1./1e3)];
px  =   [min(D.vx(1:N.nz1,:),[],2)',...
    fliplr(max(D.vx(1:N.nz1,:),[],2)')];
patch(px,py,1,'FaceColor',[.8 .8 .8],'EdgeColor','none');
alpha(0.4); hold on; box on
plot(D.vx_mean.*100.*(60*60*24*365.25),...
    M.z1'./1e3,'k','LineWidth',2)
plot(D.vx_ana.*100.*(60*60*24*365.25),...
    M.z1./1e3,'y--')
xlabel('$$\langle v_x \rangle$$ [ cm/a ]','Interpreter','latex')
ylabel('z [ km ]','Interpreter','latex')
title('Horizontal velocity','Interpreter','latex')
set(gca,'FontWeight','Bold','FontSize',15,...
    'TickLabelInterpreter','latex')
ax3=subplot(2,2,4);
plot(D.Deltav,...
    M.z1'./1e3,'k','LineWidth',2)
xlabel('$$\Delta v $$','Interpreter','latex')
ylabel('z [ km ]','Interpreter','latex')
title('Root mean square','Interpreter','latex')
set(gca,'FontWeight','Bold','FontSize',15,...
    'TickLabelInterpreter','latex')
% ======================================================================= %
%% ======================== Save final figure =========================== %
switch Pl.savefig
    case 'yes'
        saveas(figure(1),[M.ModDir,'/Channelflow'],'png')
end
% ======================================================================= %
%% ====================== Clear path structure ========================== %
if strcmp(getenv('OS'),'Windows_NT')
    rmpath('..\..\StokesProblem')
    rmpath('..\..\SetUp')
else
    rmpath('../../StokesProblem')
    rmpath('../../SetUp')
end
% ======================================================================= %
% ======================================================================= %
% =============================== END =================================== %
% ======================================================================= %


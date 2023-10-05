function [T,Py] = ContinentalGeotherm_1D(T,M,N,Py,t,plotparam)
% Function to calculate the 1D geotherm for a continental lithosphere.    %
% Temperature is calculated by solving the 1-D heat equation assuming     %
% variable thermal parameters and a radiogenic heat source.               %
% The equation is solved using a proper conserving finite difference      %
% scheme.                                                                 %
% The vertical axis is pointing downwards in negative direction.          %
%                                                                         %
% Input:                                                                  %
%   A list of structural parameters M,N,Py,t,plotparam                    %
%   M - contains:                                                         %
%       H   -   Model depth [ m ]                                         %
%       zUC -   Lower upper crust boundary [ m ]                          %
%       zLC -   Lower lower crust boundary [ m ]                          %
%   N - contains:                                                         %
%       nz  -   Number of grid points                                     %
%   Py - contains:                                                        %
%       for each layer (upper[UC] and lower[LC] crust and mantle[M])      %
%       rho -   Density [ kg/m^3 ]                                        %
%       cp  -   Specific heat [ J/kg/K ]                                  %
%       k   -   Thermal conductivity [ W/m/K ]                            %
%       H   -   Radiogentic heat source per mass [ W/kg ]                 %
%   t - contains:                                                         %
%       dtfac   -   Courant criterion                                     %
%       age     -   Lithospheric age [ Ma ]                               %
%   plotparam:                                                            %
%       0   -   no plotting                                               %
%       1   -   Plot initial and final temperature profile and thermal    %
%               parameters. The solution of the 1-D heat equation is      %
%               compared to the steady state solution (poisson equation). %
%                                                                         %
% Output:                                                                 %
%   T - contains:                                                         %
%       Tpot    -   Potential mantle temperature [ K ]                    %
%       dTadi   -   Adiabatic mantle temperature gradient [ K/km ]        %
%       T0      -   Surface temperature [ K ]                             %
%       T1      -   Bottom temperature [ K ]                              %
%       T       -   Tempeture profile [ K ]                               %
%                                                                         %
% ----------------------------------------------------------------------- %
%    LF - 17.02.2023 -                                                    %
% ======================================================================= %

if nargin==0
    %% ================================================================== %
    % Use some default values if no input parameters are defined ======== %
    % Constants --------------------------------------------------------- %
    M.H         =   -200e3;             % Hight of the model [ m ]
    M.zUC       =   -10e3;              % Depth of the upper crust [ m ]
    M.zLC       =   -35e3;              % Depth of the lower crust [ m ]
    
    N.nz        =   201;                % Number of grid points
    N.dz        =   M.H/(N.nz - 1);     % Grid resolution
    
    M.z         =   (0:N.dz:M.H)';      % Depth [ m ]
    M.zc        =   (N.dz/2:N.dz:M.H-N.dz/2)';
    
    % Mantle properties
    Py.rhoM     =   3000;               % Density [ kg/m^3 ]
    Py.cpM      =   1000;               % Heat capacity [ J/kg/K ]
    Py.kM       =   2.3;                % Conductivity [ W/m/K ]
    Py.HM       =   0;                  % Heat generation rate [W/kg]; Q = rho*H;2.3e-12
    
    % Upper crust properties
    Py.rhoUC    =   2700;               % [ kg/m^3 ]
    Py.kUC      =   3.0;                % [ W/m/K ]
    Py.HUC      =   617e-12;            % [ W/kg ]
    Py.cpUC     =   Py.cpM;
    
    % Lower crust properties
    Py.rhoLC    =   2900;               % [ kg/m^3 ]
    Py.kLC      =   2.0;                % [ W/m/K ]
    Py.HLC      =   43e-12;             % [ W/kg ]
    Py.cpLC     =   Py.cpM;
    
    M.UCind     =   M.z>=M.zUC;
    M.LCind     =   M.z>=M.zLC&M.z<M.zUC;
    M.Mind      =   M.z<=M.zLC;
    
    % Initial Condition ------------------------------------------------- %
    T.Tpot  =   1315 + 273.15;                      % Potential temperautre [ K ]
    T.dTadi =   0.5;                                % Adiabatic temperature gradient [ K/km ]
    T.T0    =   273.15;                             % Surface temperature [ K ]
    T.T1    =   T.Tpot + T.dTadi.*abs(M.H./1e3);    % Bottom temperature [ K ]
    T.T     =   T.Tpot + abs(M.z/1e3)*T.dTadi;      % Initial T-profile [ K ]
    T.T(1)  =   T.T0;    
    T.ubound    =   'flux';
    T.utbf      =   -0.0133;                % c     =   -k/q -> 50 mW/m^2
    T.lbound    =   'flux';        
    T.ltbf      =   -0.0043;                % c     =   -k/q -> 10 mW/m^2
    
    % Time stability criterion ------------------------------------------ %
    t.dtfac =   0.9;                % Courant criterion
    t.age   =   1000;                 % Lithosphere age [ Ma ]
    t.tfac  =   (60*60*24*365.25);  % Seconds per year
    t.age   =   t.age.*1e6*t.tfac;  % Age in seconds
    % =================================================================== %
    
    %% Plot Initial condition =========================================== %
    %     pos     =   [115.4,342.0,1206.4,420.0];
    %     set(figure(1),'Position',pos);
    plotparam   =   1;
    fig = figure;
    clf
    subplot(1,3,1)
    plot(T.T,M.z./1e3,'LineWidth',2)
    xlabel('T [K]'); ylabel('Depth [km]'); title('T-profile')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
    % =================================================================== %
else
    %% ================================================================== %
    % Initial Condition ------------------------------------------------- %
    T.T1    =   T.Tpot + T.dTadi.*abs(M.H./1e3);    % Bottom temperature [ K ]
    T.T     =   T.Tpot + abs(M.z/1e3)*T.dTadi;      % Initial T-profile [ K ]
    T.T(1)  =   T.T0;
    % =================================================================== %
    if plotparam
        fig = figure;
        clf
        subplot(1,3,1)
        plot(T.T,M.z./1e3,'LineWidth',2)
        xlabel('T [K]'); ylabel('Depth [km]'); title('T-profile')
        set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
    end
end
%% Setup fields ========================================================= %
Py.k            =   zeros(N.nz,1);
Py.rho          =   zeros(N.nz,1);
Py.cp           =   zeros(N.nz,1);
Py.H            =   zeros(N.nz,1);

Py.rho(:)       =   Py.rhoM;
Py.rho(M.UCind) =   Py.rhoUC;
Py.rho(M.LCind) =   Py.rhoLC;

Py.k(:)         =   Py.kM;
Py.k(M.UCind)   =   Py.kUC;
Py.k(M.LCind)   =   Py.kLC;

Py.cp(:)        =   Py.cpM;
Py.cp(M.UCind)  =   Py.cpUC;
Py.cp(M.LCind)  =   Py.cpLC;

Py.H(:)         =   Py.HM;
Py.H(M.UCind)   =   Py.HUC;
Py.H(M.LCind)   =   Py.HLC;

%% Boundary condition =================================================== %
Py.H(1)         =   0;
Py.H(end)       =   0;

%% Time stability criterion ============================================= %
Py.kappa        =   max(Py.k)/min(Py.rho)/min(Py.cp);

t.dtexp         =   N.dz^2/2/Py.kappa;  % Stability criterion for explicit
t.dt            =   t.dtfac*t.dtexp;
t.nit           =   ceil(t.age./t.dt);  % Number of iterations
t.time          =   zeros(1,t.nit);     % Time array

%% Calculate 1-D temperature profile ==================================== %
for i = 1:t.nit
    if i > 1
        t.time(i)   =   t.time(i-1) + t.dt;
    end
%     [T]     =   SolveDiff1Dimplicit_vary(N,T,Py,t);
    [T]     =   SolveDiff1Dexplicit_vary(N,T,Py,t);
end
% ----------------------------------------------------------------------- %

%% Calculate heaf flow ================================================== %
T.q         =   zeros(N.nz-1,1);
for j=1:N.nz-1
    T.q(j) = -(Py.k(j+1) + Py.k(j))/2 * ...
        (T.T(j+1) - T.T(j))/N.dz;
end
T.q(1)          =   - Py.k(1)*(T.T(2)-T.T(1))/N.dz;
T.q(N.nz-1)     =   - Py.k(N.nz)*(T.T(N.nz)-T.T(N.nz-1))/N.dz;
% ======================================================================= %

%% Plot profile if requested ============================================ %
if plotparam
    [T2]    =   Diff1D_stationary(N.nz,M.H,...
        Py.k,Py.H,Py.rho,T.T0,T.T1);
    
    if nargin == 0
        % Works only for the default parameters!
        T3      =   load('Continental_Geotherm_uc_lc_2D.txt');
    end
    
    figure(fig)
    subplot(1,3,1)
    hold on
    plot(T.T,M.z./1e3,'r-','LineWidth',2)
    plot(T2,M.z./1e3,'y--','LineWidth',2)
    if nargin == 0
%         plot(T3(:,1)+273.15,T3(:,2),'m-.','LineWidth',2)
        legend('Inital',['T_{',num2str(t.age/1e6/t.tfac),'Ma}'],...
            'T_{stationary}','T_{2D}','Location','SouthWest')
    else
        legend('Inital',['T_{',num2str(t.age/1e6/t.tfac),'Ma}'],...
            'T_{stationary}','Location','SouthWest')
    end
    xlabel(' T [K] '); ylabel(' Depth [km] '); title('T-profile')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
    subplot(1,3,2)
    plot(T.q.*1e3,M.zc./1e3,'LineWidth',2)
    xlabel('q [ mW ]'); ylabel('Depth [ km ]'); title('Heat Flux')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
    subplot(1,3,3)
    hold on
    plot(Py.k,M.z./1e3,'LineWidth',2)
    plot(Py.cp./1e3,M.z./1e3,'LineWidth',2)
    plot(Py.rho,M.z./1e3,'LineWidth',2)
    plot(Py.H.*Py.rho./1e3,M.z./1e3,'LineWidth',2)
    xlabel('k; \rho; c_{p}; Q'); ylabel('Depth [ km ]'); box on
    title('Thermal parameters')
    legend('k [W/m/K]','c_{p} [kJ/kg/K]','\rho [kg/m^3]','Q [ mW/m^3 ]',...
        'Location','SouthWest')
    set(gca,'xscale','log','FontWeight','Bold','LineWidth',2,'FontSize',15)
end
% ======================================================================= %
% keyboard
end



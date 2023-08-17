function [T,Py] = OceanicGeotherm_1D(T,M,N,Py,t,plotparam)
% Function to calculate the 1D geotherm for an oceanic lithosphere.       %
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
%   N - contains:                                                         %
%       nz  -   Number of grid points                                     %
%   Py - contains:                                                        %
%       rho -   Density [ kg/m^3 ]                                        %
%       cp  -   Specific heat [ J/kg/K ]                                  %
%       k   -   Thermal conductivity [ W/m/K ]                            %
%       H   -   Volumetric radiogenic heat source [ W/kg ]                %
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
%
% ----------------------------------------------------------------------- %
%    LF - 17.02.2023 -                                                    %
% ======================================================================= %

if nargin==0
    %% ================================================================== %
    % Use some default values if no input parameters are defined ======== %
    % Constants --------------------------------------------------------- %
    M.H         =   -200e3;             % Hight of the model [ m ]
    N.nz        =   201;                % Number of grid points
    N.dz        =   M.H/(N.nz - 1);     % Grid resolution
    M.z         =   (0:N.dz:M.H)';      % Depth [ m ]
    M.zc        =   (N.dz/2:N.dz:M.H-N.dz/2);
    Py.rhoM     =   3000;               % Density [ kg/m^3 ]
    Py.cpM      =   1000;               % Heat capacity [ J/kg/K ]
    Py.kM       =   3.0;                % Conductivity [ W/m/K ]
    Py.HM       =   0;                  % Heat generation rate [W/kg]; Q = rho*H;
    % Initial Condition ------------------------------------------------- %
    T.Tpot  =   1315 + 273.15;                      % Potential temperautre [ K ]
    T.dTadi =   0.5;                                % Adiabatic temperature gradient [ K/km ]
    T.T0    =   273.15;                             % Surface temperature [ K ]
    T.T1    =   T.Tpot + T.dTadi.*abs(M.H./1e3);    % Bottom temperature [ K ]
    T.T     =   T.Tpot + abs(M.z/1e3)*T.dTadi;      % Initial T-profile [ K ]
    T.T(1)  =   T.T0;
    T.Tini  =   T.T;
    T.ubound    =   'const';
    T.utbf      =   -0.03;          % c     =   -k/q -> 90 mW/m^2
    T.lbound    =   'const';
    T.ltbf      =   -0.0033;        % c     =   -k/q -> 10 mW/m^2
    % Time stability criterion ------------------------------------------ %
    t.dtfac =   1.0;                % Courant criterion
    t.age   =   60;                 % Lithosphere age [ Ma ]
    t.tfac  =   (60*60*24*365.25);  % Seconds per year
    t.age   =   t.age.*1e6*t.tfac;  % Age in seconds
    % =================================================================== %
    
    %% Plot Initial condition ------------------------------------------- %
    plotparam   =   1;
    fig = figure;
    clf
    subplot(1,2,1)
    plot(T.T,M.z./1e3,'LineWidth',2)
    xlabel('T [K]'); ylabel('Depth [km]'); title('T-profile')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
    % =================================================================== %
else
    %% =================================================================== %
    % Initial Condition ------------------------------------------------- %
    T.T1    =   T.Tpot + T.dTadi.*abs(M.H./1e3);    % Bottom temperature [ K ]
    T.T     =   T.Tpot + abs(M.z/1e3)*T.dTadi;      % Initial T-profile [ K ]
    T.T(1)  =   T.T0;
    T.Tini  =   T.T;
    % =================================================================== %
    if plotparam
        fig = figure;
        clf
        subplot(1,2,1)
        plot(T.T,M.z./1e3,'LineWidth',2)
        xlabel('T [K]'); ylabel('Depth [km]'); title('T-profile')
        set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
    end
end
%% Setup Fields ========================================================= %
Py.rho(:)   =   Py.rhoM;
Py.cp(:)    =   Py.cpM;
Py.k(:)     =   Py.kM;
Py.H(:)     =   Py.HM;

Py.kappa    =   max(Py.k)/min(Py.rho)/min(Py.cp); % Thermal diffusivity [ m^2/s ]
% ======================================================================= %

%% Time stability criterion ============================================= %
t.dtexp     =   N.dz^2/2/Py.kappa;  % Stability criterion for explicit
t.dt        =   t.dtfac*t.dtexp;
t.nit       =   ceil(t.age./t.dt);  % Number of iterations
t.time      =   zeros(1,t.nit);     % Time array
% ======================================================================= %

%% Calculate 1-D temperature profile ==================================== %
for i = 1:t.nit
    if i > 1
        t.time(i)   =   t.time(i-1) + t.dt;
    end
    [T]     =   SolveDiff1Dimplicit_vary(N,T,Py,t);
end
% ======================================================================= %

%% Calculate heaf flow ================================================== %
T.q         =   zeros(N.nz-1,1);
if size(Py.k,1)==1
    for j=1:N.nz-1
        T.q(j) = -Py.k * ...
            (T.T(j+1) - T.T(j))/N.dz;
    end
    T.q(1)          =   - Py.k*(T.T(2)-T.T(1))/N.dz;
    T.q(N.nz-1)     =   - Py.k*(T.T(N.nz)-T.T(N.nz-1))/N.dz;
else
    for j=1:N.nz-1
        T.q(j) = -(Py.k(j+1) + Py.k(j))/2 * ...
            (T.T(j+1) - T.T(j))/N.dz;
    end
    T.q(1)          =   - Py.k(1)*(T.T(2)-T.T(1))/N.dz;
    T.q(N.nz-1)     =   - Py.k(N.nz)*(T.T(N.nz)-T.T(N.nz-1))/N.dz;
end
% ======================================================================= %

%% Plot profile if requested ============================================ %
if plotparam
%     if nargin==0 && strcmp(T.lbound,'const') && strcmp(T.ubound,'const')
    if strcmp(T.lbound,'const') && strcmp(T.ubound,'const')
        T.Tana      =   zeros(N.nz,1);
        T.Tana      =   T.Tini + ...
            (T.T0 - T.Tpot).*erfc(-M.z./(2*sqrt(t.age*Py.kappa)));
        T.Tana(1)   =   T.T0;
    end
    figure(fig)
    subplot(1,2,1)
    hold on
    plot(T.T,M.z./1e3,'r-','LineWidth',2)
%     if nargin==0 && strcmp(T.lbound,'const') && strcmp(T.ubound,'const')
    if strcmp(T.lbound,'const') && strcmp(T.ubound,'const')
        plot(T.Tana,M.z/1e3,'--','LineWidth',2)
        legend('Initial',['T_{',num2str(t.age/1e6/t.tfac),'Ma}'],...
            'T_{HSCM}','Location','SouthWest')
    else
        legend('Initial',['T_{',num2str(t.age/1e6/t.tfac),'Ma}'],...
            'Location','SouthWest')
    end
    xlabel('T [ K ]'); ylabel('Depth [ km ]')
    subplot(1,2,2)
    plot(T.q.*1e3,M.zc./1e3,'LineWidth',2)
    xlabel('q [ mW ]'); ylabel('Depth [ km ]'); title('Heat flux')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15)
end
% ======================================================================= %

end
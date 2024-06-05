function [M,N,D,T,S] = ScaleParameters(B,M,Py,N,D,T)
%% Define Scaling Parameters ============================================ %
S.Hsc       =   abs(M.H);                   % Length scale
S.tsc       =   S.Hsc^2/Py.kappa;           % Time scale
S.vsc       =   Py.kappa/S.Hsc;             % Velocity scale
if isfield(Py,'eparam')
    S.eta       =   Py.eta0;
    S.tausc     =   S.eta*Py.kappa/S.Hsc^2;     % Stress scale
end
S.Tsc       =   Py.DeltaT;                  % Temperature scale

S.Qsc       =   Py.cp*S.Tsc*Py.kappa/S.Hsc^2;     % Heat flow scale

% Scale lenght dimensions ----------------------------------------------- %
M.H         =   M.H/S.Hsc;
M.L         =   M.L/S.Hsc;
M.X         =   M.X./S.Hsc;
M.Z         =   M.Z./S.Hsc;
M.X1        =   M.X1./S.Hsc;
M.Z1        =   M.Z1./S.Hsc;
M.x         =   M.x./S.Hsc;
M.z         =   M.z./S.Hsc;
M.x1        =   M.x1./S.Hsc;
M.z1        =   M.z1./S.Hsc;
N.dx        =   N.dx/S.Hsc;
N.dz        =   N.dz/S.Hsc;
% ----------------------------------------------------------------------- %

% Scale time dimensions ------------------------------------------------- %
if ~isempty(T)
    T.tmax      =   T.tmax/S.tsc;
end
% ----------------------------------------------------------------------- %

% Scale stokes dimensions ----------------------------------------------- %
if isfield(B,'AdvMethod')
    D.vx        =   D.vx./S.vsc;
    D.vz        =   D.vz./S.vsc;
    D.vxi       =   D.vxi./S.vsc;
    D.vzi       =   D.vzi./S.vsc;
    D.P         =   D.P./S.tausc;
    D.Pi        =   D.Pi./S.tausc;
    D.eta       =   D.eta./Py.eta0;
end
% ----------------------------------------------------------------------- %

% Scale thermal dimensions ---------------------------------------------- %
if ~isempty(D.T)
    D.T         =   (D.T-D.T(1,1))./S.Tsc;
    D.Q         =   D.Q./S.Qsc;
end
% ----------------------------------------------------------------------- %
end
% ----------------------------------------------------------------------- %

function [Py,D,ID,M,N,T,A,Pl] = SetUpFields(Py,B,N,M,T,Pl)
ID          =   [];
if ~isfield(B,'AdvMethod')
    B.AdvMethod     =   'none';
    B.Aparam        =   'none';
end

%% ------------------- Berechnung Numerischer Parameter ----------------- %
M.H         =   M.H*1000;       %   jetzt in [m]
M.L         =   -M.xmax*M.H;    %   Model Laenge  - [m]

N.nz1       =   (N.nz-1);
N.nx1       =   N.nx-1;
N.dz        =   M.H/(N.nz1);
N.dx        =   M.L/(N.nx1);

% Koordinaten des regulaeren Gitters
M.x         =   (0:N.dx:M.L);
M.z         =   (0:N.dz:M.H);
[M.X,M.Z]   =   meshgrid(M.x,M.z);

% Koordinaten des Zentrums eines FD-Elements
M.x1        =   (N.dx/2:N.dx:(M.L-N.dx/2));
M.z1        =   (N.dz/2:N.dz:(M.H-N.dz/2));
[M.X1,M.Z1] =   meshgrid(M.x1,M.z1);
A           =   [];

N.beenhere  =   0;
clear M.x M.z M.x1 M.z1
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Physikalische Parameter auf dem regularen Gitter
D.rho      =   zeros(N.nz,N.nx);
if strcmp(B.Aparam,'comp')
    D.C        =   zeros(N.nz,N.nx);
end
D.T        =    zeros(N.nz,N.nx);
if isfield(Py,'Q0')
    D.Q    =    ones(N.nz,N.nx).*Py.Q0;
else
    D.Q    =    zeros(N.nz,N.nx);
end

if isfield(Py,'tparam')
    switch Py.tparam
        case 'variable'
            D.k     = zeros(N.nz,N.nx);
    end
end

if isfield(Py,'eparam')
    % "Staggered" Gitterpunkte
    D.P         =   zeros(N.nz,N.nx);
    D.vx        =   zeros(N.nz,N.nx);
    D.vz        =   zeros(N.nz,N.nx);
    
    D.Pi        =   zeros(N.nz,N.nx);
    D.vxi       =   zeros(N.nz,N.nx);
    D.vzi       =   zeros(N.nz,N.nx);
end

if ~isempty(T)
    D.meanT     =   zeros(N.nz,T.itmax);
    D.meanV     =   zeros(T.itmax,1);
    D.meanV2    =   zeros(T.itmax,1);
    D.Nus       =   zeros(T.itmax,1);
    T.time      =   zeros(T.itmax,1);
end

if isfield(Py,'eparam')
    D.eta       =   zeros(N.nz,N.nx);
    
    % Interpolierte Parameter auf dem regulaeren Gitter
    ID.vx       =   zeros(N.nz,N.nx);
    ID.vz       =   zeros(N.nz,N.nx);
    ID.v        =   zeros(N.nz,N.nx);
    ID.vxo      =   [];
    ID.vzo      =   [];
    ID.indi     =   2:N.nz1;
    ID.indj     =   2:N.nx1;
end

% Visualization stuff --------------------------------------------------- %
Version             =   7;
Pl.TColMapName      =   'lajolla';
Pl.VColMapName      =   'imola';
Pl.EColMapName      =   'lapaz';
Pl.rhoColMapName    =   'turku';
Pl.vikColMapName    =   'vik';
Pl.tauColMapName    =   'nuuk';
Pl.epsColMapName    =   'batlowW';
Pl.PColMapName      =   'hawaii';

if strcmp(getenv('OS'),'Windows_NT')
    Pl.SciColDir    =   ['D:\Users\lukas\Numerics\BACKUP\',...
        'progs\src\MATLAB\FDCSGm\',...
        'ScientificColourMaps',num2str(Version)];
    Pl.Slash        =   '\';
else
    Pl.SciColDir    =   ['/home/external_homes/lufuchs/',...
        'progs/src/MATLAB/FDCSGm/',...
        '/ScientificColourMaps',num2str(Version)];
    Pl.Slash        =   '/';
end

Pl.TColMap      =   [Pl.SciColDir,Pl.Slash,Pl.TColMapName,Pl.Slash,...
    Pl.TColMapName,'.mat'];
Pl.VColMap      =   [Pl.SciColDir,Pl.Slash,Pl.VColMapName,Pl.Slash,...
    Pl.VColMapName,'.mat'];
Pl.EColMap      =   [Pl.SciColDir,Pl.Slash,Pl.EColMapName,Pl.Slash,...
    Pl.EColMapName,'.mat'];
Pl.rhoColMap    =   [Pl.SciColDir,Pl.Slash,Pl.rhoColMapName,Pl.Slash,...
    Pl.rhoColMapName,'.mat'];
Pl.vikColMap    =   [Pl.SciColDir,Pl.Slash,Pl.vikColMapName,Pl.Slash,...
    Pl.vikColMapName,'.mat'];
Pl.tauColMap    =   [Pl.SciColDir,Pl.Slash,Pl.tauColMapName,Pl.Slash,...
    Pl.tauColMapName,'.mat'];
Pl.epsColMap    =   [Pl.SciColDir,Pl.Slash,Pl.epsColMapName,Pl.Slash,...
    Pl.epsColMapName,'.mat'];
Pl.PColMap      =   [Pl.SciColDir,Pl.Slash,Pl.PColMapName,Pl.Slash,...
    Pl.PColMapName,'.mat'];

load(Pl.TColMap); Pl.lajolla    = lajolla;
load(Pl.VColMap); Pl.imola      = imola;
load(Pl.EColMap); Pl.lapaz      = lapaz;
load(Pl.rhoColMap); Pl.oslo     = turku;
load(Pl.vikColMap); Pl.vik      = vik;
load(Pl.tauColMap); Pl.nuuk     = nuuk;
load(Pl.epsColMap); Pl.batlowW  = batlowW;
load(Pl.PColMap); Pl.hawaii  = hawaii;

end

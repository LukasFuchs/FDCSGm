function [Anew] = SemiLagAdvection2DSc(D,M,A,dt)
% Funktion zur Loesung der Advektionsgleichung
%
%               dA/dt = -vx * ( dA/dx ) - vz * ( dA/dz )
%
% in einer zweidimensionalen Umgebung mit Hilfe des Semi-Lagrangian Schema,
% und einer Zentralen-Punkt Iteration, d.h.
%
% ======================================================================= %

% mid-point iteration scheme -------------------------------------------- %
if isempty(D.vxo)
    % Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
    % aktuelle Geschwindigkeit auf die alte Geschwindigkeit
    % uebertragen.
    D.vxo   =   D.vx;
    D.vzo   =   D.vz;
end

% Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --------------- %
vxm     = 0.5.*(D.vxo + D.vx);
vzm     = 0.5.*(D.vzo + D.vz);

% Initialisierung der Geschwindigkeit fuer die Iteration ---------------- %
vxi     = D.vx;
vzi     = D.vz;

% Iteration ------------------------------------------------------------- %
for k = 1:10
    xp  = M.X - 0.5*dt.*vxi;
    zp  = M.Z - 0.5*dt.*vzi;
    
    vxi = interp2(M.X,M.Z,vxm,xp,zp,'linear');
    vzi = interp2(M.X,M.Z,vzm,xp,zp,'linear');
    
    vxi(isnan(vxi)) = vxm(isnan(vxi));
    vzi(isnan(vzi)) = vzm(isnan(vzi));
end
xp      =   M.X - dt.*vxi;
zp      =   M.Z - dt.*vzi;

Anew    =   interp2(M.X,M.Z,A,xp,zp,'cubic');

Anew(isnan(Anew)) = A(isnan(Anew));

end
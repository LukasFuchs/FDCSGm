function [D,B,Ma,Py] = SetUpInitialConditions(D,Py,M,N,B)

Ma  = [];

%% Randbedingungen ------------------------------------------------------ %
B.B0c       =   N.nz*2+2;
D.Pi(B.B0c) =   0;
% ----------------------------------------------------------------------- %

%% Temperaturanfangsbedingungen ----------------------------------------- %
switch B.Tini
    case 'gaussian'
        % Gaussche Temperatur Anomalie ---------------------------------- %
        Ampl    =   10;            % Amplitude der Anomalie [ K ]
        sigma   =   M.L/10;         % Breite der Anomalie
        T0      =   1000;           % Hintergrund Temperatur [ K ]
        
        D.T     =   T0 + ...
            Ampl.*exp(-((M.X-0.25*M.L).^2+(M.Z-0.5*M.H).^2)./sigma^2);       
        
    case 'block'
        % Hintergrund Temperatur ---------------------------------------- %
        Tb      =   1000;           % [ K ]
        
        % Bereich der Temperatur Anomalie ------------------------------- %
        xTl     =   M.L/8;
        xTr     =   xTl + M.L/5;
        zTu     =   M.H/2 + M.H/10;
        zTo     =   M.H/2 - M.H/10;
        
        Ta      =   1010;           % [ K ]       
        
        % Anfangstemperatur Verteilung ---------------------------------- %
        D.T    =   Tb.*ones(N.nz,N.nx);
        D.T(M.z>=zTu&M.z<=zTo,M.x>=xTl&M.x<=xTr) = Ta;
    case 'const'
        D.T    =   ones(N.nz,N.nx).*1000;           % [ K ]
    case 'linear'
        Ttop        =   273;
        Tgrad       =   0.5;                            % [ K/km ]
        D.T         =   -Tgrad*M.Z./1e3 + Ttop;
        B.bhf       =   D.T(N.nz,1); 
        B.thf       =   D.T(1,1); 
        Py.DeltaT   =   B.bhf - B.thf; 
    case 'linano'
        Ttop        =   273;
        Tgrad       =   0.5;                            % [ K/km ]
        D.T         =   -Tgrad*M.Z./1e3 + Ttop;
        B.bhf       =   D.T(N.nz,1); 
        B.thf       =   D.T(1,1); 
        Py.DeltaT   =   B.bhf - B.thf; 
        
        % Gaussche Temperatur Anomalie ---------------------------------- %
        % Amplitude der Anomalie [ K ]
        Ampl    =   D.T(find(M.Z(:,1)==0.5*M.H),1)*0.1; 
        sigma   =   M.L/10;         % Breite der Anomalies
        zano    =   0.5.*M.H;       % Tiefe des Zentrums der Anomalie
        
        T0      =   Ampl.*...
            exp(-((M.X-0.25*M.L).^2+(M.Z-zano).^2)./sigma^2);
        
        D.T     = D.T + T0; 
        
    otherwise
        error('Initial Temperaturefield not defined!')
end

switch B.btbc
    case 'const'
        D.T(N.nz,:) = B.bhf;                        % [ K ]
end
switch B.ttbc
    case 'const'
        D.T(1,:)    = B.thf;                        % [ K ]
end
% ----------------------------------------------------------------------- %

%% Anfangsgeschwindigkeitsfeld ========================================== %
switch B.IniFlow
    case 'RigidBody'
        % Starre Rotation
        D.vx    =   -B.FlowFac.*((M.Z-N.dz/2)-M.H/2)./M.H;       % [ cm/a ]
        D.vz    =   B.FlowFac.*((M.X-N.dx/2)-M.L/2)./M.H;      % [ cm/a ]
        
        D.vx    =   D.vx/(100*(60*60*24*365.25));   % [ m/s ]
        D.vz    =   D.vz/(100*(60*60*24*365.25));   % [ m/s ]
        
    case 'ShearCell'
        % Zelle mit einfacher Scherdehnung
        D.vx    =   B.FlowFac.*(sin(pi.*M.X./M.L).*cos(-pi.*(M.Z-N.dz/2)./M.H));
        D.vz    =   -B.FlowFac.*(cos(pi.*(M.X-N.dx/2)./M.L).*sin(-pi.*M.Z./M.H));
        
        D.vx    =   D.vx/(100*(60*60*24*365.25));   % [ m/s ]
        D.vz    =   D.vz/(100*(60*60*24*365.25));   % [ m/s ]
    case 'none'
        D.vx    =   zeros(N.nz,N.nx);
        D.vz    =   zeros(N.nz,N.nx);
    case 'constPlate'
%         D.vzi(1,1:N.nx1)    = -1;        % cm/a
%         D.vzi(1:N.nz1,N.nx) = 1;        % cm/a        
        
%         D.vx   =   D.vx./(100*(60*60*24*365.25));   % [ m/s ]
%         D.vzi   =   D.vzi./(100*(60*60*24*365.25));   % [ m/s ]
    otherwise
        error('Flowfield not defined!')
end
switch B.AdvMethod
    case 'none'
        D.vx    = zeros(N.nz,N.nx);
        D.vz    = zeros(N.nz,N.nx);
end
% ----------------------------------------------------------------------- %

%% Anfangsdichte -------------------------------------------------------- %
% Zustandsgleichung 
D.rho       =   Py.rho0.*(1-Py.alpha.*(D.T-D.T(1,1)));
% ----------------------------------------------------------------------- %

%% Anfangsviskositaet --------------------------------------------------- %
switch Py.eparam
    case 'variable'                
        D.eta       =   ones(N.nz,N.nx); 
        
        D.eta       =   Py.eta0.*...
            exp( -Py.b.*((D.T-D.T(1,1))./(D.T(nz,1)-D.T(1,1))) + ...
            Py.c.*M.Z./M.H);
end
% ----------------------------------------------------------------------- %

end
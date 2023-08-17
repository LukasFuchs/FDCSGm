function [T,D,B,M,Ma,Py] = SetUpInitialConditions(T,D,Py,M,N,B)

Ma  = [];

%% Randbedingungen ------------------------------------------------------ %
B.B0c       =   N.nz*2+2;
if isfield(Py,'eparam')
    D.Pi(B.B0c) =   0;
else
    B.EtaIni    =   'none';
end
% ----------------------------------------------------------------------- %

%% Maximale Zeit -------------------------------------------------------- %
if ~isempty(T)
    T.tmax      =   T.tmaxini.*1e6*(365.25*24*60*60);      % in [ s ]
end
% ----------------------------------------------------------------------- %

%% Temperaturanfangsbedingungen ----------------------------------------- %
if isfield(B,'Tini')
    switch B.Tini
        case 'none'
            D.T     = [];
        case 'ellipse'
            Tb      = Py.DeltaT/2;
            Ta      = Tb + Py.DeltaT*0.2;
            % Bereich der Temperatur Anomalie --------------------------- %
            % Halbachsen der Ellipse in m
            a       = -M.H/2;
            b       = -M.H/6;
            
            % Mittelpunkt der Ellipse
            xc      = M.L/4;
            zc      = M.H/2 - b;
            
            Elleq   = (M.X-xc).^2./a^2 + (M.Z-zc).^2./b^2;
            
            ind     = Elleq < 1;
            
            D.T(ind)    = Ta;
            D.T(~ind)   = Tb;
            
        case 'gaussian'
            % Gaussche Temperatur Anomalie ------------------------------ %
            B.Tsigma    =   (B.Tsigma/100)*M.L; % Breite der Anomalie [ m ]
            
            D.T         =   B.T0 + ...
                B.TAmpl .*exp( -( (M.X-0.5*M.L).^2 + (M.Z-0.5*M.H).^2 ) ./...
                ( 2*B.Tsigma^2/pi ));
            
            D.Tana      =   D.T;
            
            D.TProfile  =   zeros(N.nz,T.itmax);
            D.Tmax      =   zeros(1,T.itmax);
            D.TProfilea =   zeros(N.nz,T.itmax);
            D.Tmaxa     =   zeros(1,T.itmax);
            
            %                 D.T     =   B.T0 + ...
            %             B.TAmpl.*exp(-((M.X-0.25*M.L).^2+(M.Z-0.5*M.H).^2)./sigma^2);
            
        case 'block'
            % Hintergrund Temperatur ------------------------------------ %
            Tb      =   1000;           % [ K ]
            
            % Bereich der Temperatur Anomalie --------------------------- %
            xTl     =   M.L/8;
            xTr     =   xTl + M.L/5;
            zTu     =   M.H/2 + M.H/10;
            zTo     =   M.H/2 - M.H/10;
            
            Ta      =   1020;           % [ K ]
            
            % Anfangstemperatur Verteilung ------------------------------ %
            D.T    =   Tb.*ones(N.nz,N.nx);
            D.T(M.z>=zTu&M.z<=zTo,M.x>=xTl&M.x<=xTr) = Ta;
        case 'const'
            D.T    =   ones(N.nz,N.nx).*1000;               % [ K ]
        case 'linear'
            Ttop        =   273;
            Tgrad       =   0.5;                            % [ K/km ]
            D.T         =   -Tgrad*M.Z./1e3 + Ttop;
            B.bhf       =   D.T(N.nz,1);
            B.thf       =   D.T(1,1);
            Py.DeltaT   =   B.bhf - B.thf;
        case 'ContGeotherm'
            
        case 'linano'
            Ttop        =   273;
            Tgrad       =   0.5;                            % [ K/km ]
            D.T         =   -Tgrad*M.Z./1e3 + Ttop;
            B.bhf       =   D.T(N.nz,1);
            B.thf       =   D.T(1,1);
            Py.DeltaT   =   B.bhf - B.thf;
            
            % Gaussche Temperatur Anomalie ------------------------------ %
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
end

if isfield(Py,'tparam')
    if ~strcmp(Py.tparam,'none')
        switch B.btbc
            case 'const'
                D.T(N.nz,:)         = B.bhf;                % [ K ]
        end
        switch B.ttbc
            case 'const'
                D.T(1,:)            = B.thf;                % [ K ]
        end
        switch B.ltbc
            case 'const'
                D.T(2:N.nz1,1)      = B.lhf;                % [ K ]
        end
        switch B.rtbc
            case 'const'
                D.T(2:N.nz1,N.nx)   = B.rhf;                % [ K ]
        end
    else
        B.Tini 	= 	'none';
    end
end
% ======================================================================= %
%% Anfangsgeschwindigkeitsfeld ========================================== %
if isfield(B,'IniFlow')
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
            
        case 'PureShear'
            D.vxi(1:N.nz1,1)        = (M.X(1:N.nz1,end)-M.L/2).*B.ebg;    % left
            D.vxi(1:N.nz1,N.nx)     = -(M.X(1:N.nz1,end)-M.L/2).*B.ebg;   % right
            D.vxi(1,2:N.nx1)        = -(M.X(1,2:N.nx1)-M.L/2).*B.ebg;      % top
            D.vxi(N.nz1,2:N.nx1)    = -(M.X(N.nz1,2:N.nx1)-M.L/2).*B.ebg;  % bottom
            
            D.vzi(2:N.nz1,1)        = (M.Z(2:N.nz1,1)-M.H/2).*B.ebg;      % left
            D.vzi(2:N.nz1,N.nx1)    = (M.Z(2:N.nz1,N.nx1)-M.H/2).*B.ebg;  % right
            D.vzi(1,1:N.nx1)        = (M.Z(1,1:N.nx1)-M.H/2).*B.ebg;      % top
            D.vzi(N.nz,1:N.nx1)     = (M.Z(N.nz,1:N.nx1)-M.H/2).*B.ebg;   % bottom
            
            %         D.vxi   =   D.vxi./(100*(60*60*24*365.25));   % [ m/s ]
            %         D.vzi   =   D.vzi./(100*(60*60*24*365.25));   % [ m/s ]
        case 'SimpleShear'
            D.vxi(1:N.nz1,1)        = (M.Z1(1:N.nz1,end)-M.H/2).*B.ebg;    % left
            D.vxi(1:N.nz1,N.nx)     = (M.Z1(1:N.nz1,end)-M.H/2).*B.ebg;    % right
            D.vxi(1,2:N.nx1)        = (M.Z1(1,2:N.nx1)-M.H/2).*B.ebg;      % top
            D.vxi(N.nz1,2:N.nx1)    = (M.Z1(N.nz1,2:N.nx1)-M.H/2).*B.ebg;  % bottom
            
            D.vzi(2:N.nz1,1)        = 0;    % left
            D.vzi(2:N.nz1,N.nx1)    = 0;    % right
            D.vzi(1,1:N.nx1)        = 0;    % top
            D.vzi(N.nz,1:N.nx1)     = 0;    % bottom
        case 'constPlate'
            error('Not defined yet!')
        case 'none'
            D.vxi    =   zeros(N.nz,N.nx);
            D.vzi    =   zeros(N.nz,N.nx);
        otherwise
            error('Flowfield not defined!')
    end
end
% ----------------------------------------------------------------------- %

%% Anfangsdichte -------------------------------------------------------- %
% Zustandsgleichung
if isfield(Py,'tparam')
    if strcmp(Py.tparam,'variable')
        D.rho       =   Py.rho0.*(1-Py.alpha.*(D.T-D.T(1,1)));
    else
        D.rho       =   Py.rho0.*ones(N.nz,N.nx);
    end
end
% ----------------------------------------------------------------------- %

%% Anfangsviskositaet --------------------------------------------------- %
switch B.EtaIni
    case 'none'
        D.eta 	    =   [];
    case 'tdep'
        D.eta       =   Py.eta0.*...
            exp( -Py.b.*((D.T-D.T(1,1))./(D.T(N.nz,1)-D.T(1,1))) + ...
            Py.c.*M.Z./M.H);
    case 'block'
        % Block geometry
        xL          =   2/5*M.L;
        xR          =   3/5*M.L;
        zT          =   0.1*M.H;
        zB          =   0.3*M.H;
        M.ind       =   M.X>xL&M.X<xR&M.Z>zB&M.Z<zT;
        
        % Tracer Initialisierung
        nmxx        =   (N.nx-1)*N.nmx;
        nmzz        =   (N.nz-1)*N.nmz;
        dmx         =   M.L/(nmxx-1);
        dmz         =   M.H/(nmzz-1);
        xm          =   linspace(0,M.L-dmx,nmxx);
        zm          =   linspace(M.H-dmz,0,nmzz);
        
        [XM,ZM]     =   meshgrid(xm,zm);
        XM          =   XM + rand(nmzz,nmxx)*dmx;
        Ma.XM       =   reshape(XM,[nmzz*nmxx,1]);
        ZM          =   ZM + rand(nmzz,nmxx)*dmz;
        Ma.ZM       =   reshape(ZM,[nmzz*nmxx,1]);
        clear XM ZM
        Ma.C        =   ones(nmzz*nmxx,1);
        
        ind         = Ma.XM>xL & Ma.XM<xR & Ma.ZM > zB & Ma.ZM < zT;
        
        Ma.C(ind) = 2;
        
        Ma.rho      =   [Py.rho0,Py.rho1];
        Ma.eta      =   [Py.eta0,Py.eta1];
        
        % Interpolate from tracers to grid ------------------------------ %
        [~,D.eta]   =   TracerInterp(Ma,D.eta,[],M.X,M.Z,'from','eta');
        [~,D.rho]   =   TracerInterp(Ma,D.rho,[],M.X,M.Z,'from','rho');
        
        if ~isfield(B,'AdvMethod')
            Ma  =   [];
        else            
            switch lower(B.AdvMethod)
                case 'semi-lag'
                    Ma  = [];
            end
        end
        
    case 'ellipse'
        % Tracer Initialisierung
        nmxx        =   (N.nx-1)*N.nmx;
        nmzz        =   (N.nz-1)*N.nmz;
        dmx         =   M.L/(nmxx-1);
        dmz         =   M.H/(nmzz-1);
        xm          =   linspace(0,M.L-dmx,nmxx);
        zm          =   linspace(M.H-dmz,0,nmzz);
        
        [XM,ZM]     =   meshgrid(xm,zm);
        XM          =   XM + rand(nmzz,nmxx)*dmx;
        Ma.XM       =   reshape(XM,[nmzz*nmxx,1]);
        ZM          =   ZM + rand(nmzz,nmxx)*dmz;
        Ma.ZM       =   reshape(ZM,[nmzz*nmxx,1]);
        
        Ma.C        =   ones(nmzz*nmxx,1);
        
        % Bereich der Inclusion --------------------------------- %
        % Halbachsen der Ellipse in m
        a       = B.EllA;
        b       = B.EllB;
        
        alpha   = B.RotAng;
        
        % Mittelpunkt der Ellipse
        xc      = M.L/2;
        zc      = M.H/2;
        
        x_ell   = (Ma.XM-xc).*cosd(alpha) + (Ma.ZM-zc).*sind(alpha);
        z_ell   = -(Ma.XM-xc).*sind(alpha) + (Ma.ZM-zc).*cosd(alpha);
        
        Elleq   = (x_ell./a).^2 + (z_ell./b).^2;
        
        ind     =   Elleq < 1;
        
        Ma.C(ind) = 2;
        
        % Ma.rho      =   [Py.rho0,Py.rho1];
        Ma.eta      =   [Py.eta0,Py.eta1];
        
        % Interpolate from tracers to grid ------------------------------ %
        [~,D.eta]  = TracerInterp(Ma,D.eta,[],M.X,M.Z,'from','eta');
        %         [~,Py.rho]  = TracerInterp(Ma,Py.eta,[],M.X,M.Z,'from','rho');
        %                 D.eta(ind)  =   Py.eta1;
        %                 D.eta(~ind) =   Py.eta0;
        Ma  = [];
    case 'RTI'
        % Tracer Initialisierung
        nmxx        =   (N.nx-1)*N.nmx;
        nmzz        =   (N.nz-1)*N.nmz;
        dmx         =   M.L/(nmxx-1);
        dmz         =   M.H/(nmzz-1);
        xm          =   linspace(0,M.L-dmx,nmxx);
        zm          =   linspace(M.H-dmz,0,nmzz);
        
        [XM,ZM]     =   meshgrid(xm,zm);
        XM          =   XM + rand(nmzz,nmxx)*dmx;
        ZM          =   ZM + rand(nmzz,nmxx)*dmz;
        
        C           =   zeros(nmzz,nmxx);
        
        % Layer interface
        deltaAm     =   cos(2*pi*((xm - 0.5*M.L)/(B.lambda*1e3)))*B.deltaA;
        
        for j = 1:nmzz
            C(j,:)  =   ZM(j,:)  <  M.H/2 + deltaAm;
        end
        C           =   C + 1;
        
        Ma.ZM       =   reshape(ZM,[nmzz*nmxx,1]);
        Ma.XM       =   reshape(XM,[nmzz*nmxx,1]);
        Ma.C        =   reshape(C,[nmzz*nmxx,1]);
        
        Ma.c        =   [1,2];
        
        switch lower(B.AdvMethod)
            case 'semi-lag'
                % Interpolate from tracers to grid ---------------------- %
                [~,D.C]     =   TracerInterp(Ma,D.C,[],M.X,M.Z,'from','comp');
                D.C         =   round(D.C);
                Ma  = [];
                
                D.rho(D.C<=1)   =    Py.rho0;
                D.rho(D.C==2)   =    Py.rho1;
                D.eta(D.C<=1)   =    Py.eta0;
                D.eta(D.C==2)   =    Py.eta1;
            case 'tracers'
                Ma.rho      =   [Py.rho0 Py.rho1];
                Ma.eta      =   [Py.eta0 Py.eta1];
                % Interpolate from the tracers to the grid -------------- %
                [~,D.rho]   =   TracerInterp(Ma,D.rho,[],M.X,M.Z,'from','rho');
                [~,D.eta]   =   TracerInterp(Ma,D.eta,[],M.X,M.Z,'from','eta');
        end
end
% ----------------------------------------------------------------------- %
end

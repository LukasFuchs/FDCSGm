function [D,Ma,ID] = Advection(it,N,B,D,ID,Py,dt,M,Ma)
% ======================================================================= %
% Generell function to call the different advection schemes.              %
% ======================================================================= %
switch lower(B.AdvMethod)   
    case 'upwind'
                D.T             =   UpwindAdvection2D(ID.vx,ID.vz,D.T,...
                                        N.dx,N.dz,dt);
    case 'slf'
                D.Tnew          =   SLFAdvection2D(ID.vx,ID.vz,...  
                                        D.Told,D.T,N.dx,N.dz,dt);                                
                D.Told          =   D.T; 
                D.T             =   D.Tnew;
    case 'semi-lag'        
        switch lower(B.Aparam)
            case 'rho'
                if (it==1)
                    [D.rho]     =   SemiLagAdvection2D(ID,M,D.rho,dt);                    
                    % Speicher die alte Geschwindigkeit
                    ID.vxo      =   ID.vx;
                    ID.vzo      =   ID.vz;
                else
                    [D.rho]     =   SemiLagAdvection2D(ID,M,D.rho,dt);
                end
            case 'temp'
                if (it==1)
                    [D.T]       =   SemiLagAdvection2D(ID,M,D.T,dt);                    
                    % Speicher die alte Geschwindigkeit
                    ID.vxo      =   ID.vx;
                    ID.vzo      =   ID.vz;
                else
                    [D.T]       =   SemiLagAdvection2D(ID,M,D.T,dt);
                end
            case 'comp'
                if (it==1)
                    [D.C]       =   SemiLagAdvection2D(ID,M,D.C,dt);                    
                    % Speicher die alte Geschwindigkeit
                    ID.vxo      =   ID.vx;
                    ID.vzo      =   ID.vz;
                else
                    [D.C]       =   SemiLagAdvection2D(ID,M,D.C,dt);
                end
                % Update viscosity and density ========================== %
                if (mod(it,25)==0)
                    D.C         =   round(D.C);
                end
                D.rho(D.C<=1.5) =    Py.rho0;   D.rho(D.C>=1.5) =    Py.rho1;
                D.eta(D.C<=1.5) =    Py.eta0;   D.eta(D.C>=1.5) =    Py.eta1;
        end
    case 'tracers'
        % Advect tracers with Runge-Kutta 4th order --------------------- %
        [Ma.XM,Ma.ZM] = ...
            AdvectMarker2D(M,Ma.XM,Ma.ZM,dt,ID.vx,ID.vz);
        switch lower(B.Aparam)
            case 'temp'
                % Interpolate from the tracers to the grid -------------- %
                [~,D.T]     =   TracerInterp(Ma,D.T,[],M.X,M.Z,'from','temp');
            case 'comp'
                % Interpolate from the tracers to the grid -------------- %
                [~,D.rho]   =   TracerInterp(Ma,D.rho,[],M.X,M.Z,'from','rho');
                [~,D.eta]   =   TracerInterp(Ma,D.eta,[],M.X,M.Z,'from','eta');
        end
    otherwise 
        error('Advection method unknown!')
end
end
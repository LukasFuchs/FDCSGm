function [Ma,PG] = TracerInterp(Ma,PG,PGo,X,Z,direction,param)

nmr     =   size(Ma.XM,1);

nz      =   size(PG,1);
nx      =   size(PG,2);

dx      =   max(diff(X(1,:)));
dz      =   max(diff(Z(:,1)));

j       =   floor(Ma.XM./dx) + 1;
i       =   floor(Ma.ZM./dz) + 1;

switch param
    case 'comp'
        MAPROP  = Ma.c;
    case 'eta'
        MAPROP  = Ma.eta;
    case 'rho'
        MAPROP  = Ma.rho;
    case 'temp'
        MAPROP  = Ma.T;
end

switch direction
    case 'to'
        if isempty(PGo)
            % Interpolate the absolute value onto the tracers
            for k = 1:nmr
                % Normalized distance to next grid point
                dxm     = (Ma.XM(k) - X(i(k),j(k)))/dx;
                dzm     = (Ma.ZM(k) - Z(i(k),j(k)))/dz;
                
                MAPROP(k)   = PG(i(k),j(k))*(1-dxm)*(1-dzm) + ...
                    PG(i(k),j(k)+1)*dxm*(1-dzm) + ...
                    PG(i(k)+1,j(k))*dzm*(1-dxm) + ...
                    PG(i(k)+1,j(k)+1)*dxm*dzm;
            end
        else
            %             % Interpolate only the increment
            %             % Not working yet!
            %             for k = 1:nmr
            %                 % Normalized distance to next grid point
            %                 dzm     = (ZM(k) - Z(i(k),j(k)))/dz;
            %                 dxm     = (XM(k) - X(i(k),j(k)))/dx;
            %
            %                 PM(k)   = PM(k) + ...
            %                     (PG(i(k),j(k))     - PGo(i(k),j(k)))*(1-dxm)*(1-dzm) + ...
            %                     (PG(i(k),j(k)+1) - PGo(i(k),j(k)+1))*dxm*(1-dzm) + ...
            %                     (PG(i(k)+1,j(k)) - PGo(i(k)+1,j(k)))*dzm*(1-dxm) + ...
            %                     (PG(i(k)+1,j(k)+1) - PGo(i(k)+1,j(k)+1))*dxm*dzm;
            %             end
        end
        
        switch param
            case 'comp'
                Ma.c    =   MAPROP;
            case 'eta'
                Ma.eta  =   MAPROP;
            case 'rho'
                Ma.rho  =   MAPROP;
            case 'temp'
                Ma.T    =   MAPROP;
        end
        
    case 'from'
        PG  = zeros(nz,nx);
        wt  = zeros(nz,nx);
        
        switch param
            case 'temp'
                for k = 1:nmr
                    dzm     = (Ma.ZM(k) - Z(i(k),j(k)))/dz;
                    dxm     = (Ma.XM(k) - X(i(k),j(k)))/dx;
                    
                    % Upper-Left node
                    PG(i(k),j(k))     = PG(i(k),j(k)) + (1.0-dxm)*(1.0-dzm)*MAPROP(k);
                    wt(i(k),j(k))     = wt(i(k),j(k)) + (1.0-dxm)*(1.0-dzm);
                    % Lower-Left node
                    PG(i(k)+1,j(k))   = PG(i(k)+1,j(k)) + (1.0-dxm)*dzm*MAPROP(k);
                    wt(i(k)+1,j(k))   = wt(i(k)+1,j(k)) + (1.0-dxm)*dzm;
                    % Upper-Right node
                    PG(i(k),j(k)+1)   = PG(i(k),j(k)+1) + dxm*(1.0-dzm)*MAPROP(k);
                    wt(i(k),j(k)+1)   = wt(i(k),j(k)+1) + dxm*(1.0-dzm);
                    % Lower-Right node
                    PG(i(k)+1,j(k)+1) = PG(i(k)+1,j(k)+1) + dxm*dzm*MAPROP(k);
                    wt(i(k)+1,j(k)+1) = wt(i(k)+1,j(k)+1) + dxm*dzm;
                end
            case 'rho'
                for k = 1:nmr
                    dxm     = (Ma.XM(k) - X(i(k),j(k)))/dx;
                    dzm     = (Ma.ZM(k) - Z(i(k),j(k)))/dz;
                    
                    % Upper-Left node
                    PG(i(k),j(k))     = PG(i(k),j(k)) + (1.0-dxm)*(1.0-dzm)*MAPROP(Ma.C(k));
                    wt(i(k),j(k))     = wt(i(k),j(k)) + (1.0-dxm)*(1.0-dzm);
                    % Lower-Left node
                    PG(i(k)+1,j(k))   = PG(i(k)+1,j(k)) + (1.0-dxm)*dzm*MAPROP(Ma.C(k));
                    wt(i(k)+1,j(k))   = wt(i(k)+1,j(k)) + (1.0-dxm)*dzm;
                    % Upper-Right node
                    PG(i(k),j(k)+1)   = PG(i(k),j(k)+1) + dxm*(1.0-dzm)*MAPROP(Ma.C(k));
                    wt(i(k),j(k)+1)   = wt(i(k),j(k)+1) + dxm*(1.0-dzm);
                    % Lower-Right node
                    PG(i(k)+1,j(k)+1) = PG(i(k)+1,j(k)+1) + dxm*dzm*MAPROP(Ma.C(k));
                    wt(i(k)+1,j(k)+1) = wt(i(k)+1,j(k)+1) + dxm*dzm;
                end
            otherwise
                for k = 1:nmr
                    dzm     =   (Ma.ZM(k) - Z(i(k),j(k)))/dz;
                    dxm     =   (Ma.XM(k) - X(i(k),j(k)))/dx;
                    
                    % Upper-Left node
                    if(dxm <= 0.5 && dzm <= 0.5)
                        PG(i(k),j(k))     = PG(i(k),j(k)) + (1.0-dxm)*(1.0-dzm)*MAPROP(Ma.C(k));
                        wt(i(k),j(k))     = wt(i(k),j(k)) + (1.0-dxm)*(1.0-dzm);
                    end
                    % Lower-Left node
                    if(dxm<=0.5 && dzm >= 0.5)
                        PG(i(k)+1,j(k))   = PG(i(k)+1,j(k)) + (1.0-dxm)*dzm*MAPROP(Ma.C(k));
                        wt(i(k)+1,j(k))   = wt(i(k)+1,j(k)) + (1.0-dxm)*dzm;
                    end
                    % Upper-Right node
                    if(dxm >= 0.5 && dzm <= 0.5)
                        PG(i(k),j(k)+1)   = PG(i(k),j(k)+1) + dxm*(1.0-dzm)*MAPROP(Ma.C(k));
                        wt(i(k),j(k)+1)   = wt(i(k),j(k)+1) + dxm*(1.0-dzm);
                    end
                    % Lower-Right node
                    if (dxm >= 0.5 && dzm >= 0.5)
                        PG(i(k)+1,j(k)+1) = PG(i(k)+1,j(k)+1) + dxm*dzm*MAPROP(Ma.C(k));
                        wt(i(k)+1,j(k)+1) = wt(i(k)+1,j(k)+1) + dxm*dzm;
                    end
                    
                    %                     dxm     = (Ma.XM(k) - X(i(k),j(k)))/dx;
                    %                     dzm     = (Ma.ZM(k) - Z(i(k),j(k)))/dz;
                    %
                    %                     % Upper-Left node
                    %                     PG(i(k),j(k))     = PG(i(k),j(k)) + (1.0-dxm)*(1.0-dzm)*MAPROP(Ma.C(k));
                    %                     wt(i(k),j(k))     = wt(i(k),j(k)) + (1.0-dxm)*(1.0-dzm);
                    %                     % Lower-Left node
                    %                     PG(i(k)+1,j(k))   = PG(i(k)+1,j(k)) + (1.0-dxm)*dzm*MAPROP(Ma.C(k));
                    %                     wt(i(k)+1,j(k))   = wt(i(k)+1,j(k)) + (1.0-dxm)*dzm;
                    %                     % Upper-Right node
                    %                     PG(i(k),j(k)+1)   = PG(i(k),j(k)+1) + dxm*(1.0-dzm)*MAPROP(Ma.C(k));
                    %                     wt(i(k),j(k)+1)   = wt(i(k),j(k)+1) + dxm*(1.0-dzm);
                    %                     % Lower-Right node
                    %                     PG(i(k)+1,j(k)+1) = PG(i(k)+1,j(k)+1) + dxm*dzm*MAPROP(Ma.C(k));
                    %                     wt(i(k)+1,j(k)+1) = wt(i(k)+1,j(k)+1) + dxm*dzm;
                end
        end
        
        PG  = PG./wt;
        
        PG(isinf(PG))   = min(min(PG));
        
end
end

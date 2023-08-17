function [T] = SolveDiff1Dimplicit_vary(N,T,Py,t)

T0      =   T.T;

if size(Py.k,1) == 1
    k   =   Py.k.*ones(N.nz,1);
    rho =   Py.rho.*ones(N.nz,1);
    cp  =   Py.cp.*ones(N.nz,1);
else
    k   =   Py.k;
    rho =   Py.rho;
    cp  =   Py.cp;
end
if size(Py.H,1) == 1
    H   =   Py.H.*ones(N.nz,1);     % [Q] = W/m^3; [Q] = [rho*H]
else
    H   =   Py.H;
end

ind         = 2:N.nz-1;

diag        = zeros(N.nz,3);

diag(ind-1,1)   = -((k(ind)+k(ind-1))./2).*...
    t.dt/N.dz^2./rho(ind)./cp(ind);

diag(ind,2)     = (1 + t.dt/N.dz^2./rho(ind)./cp(ind).*...
    ( (k(ind+1)+k(ind))./2 + (k(ind)+k(ind-1))./2 ));

diag(ind+1,3)   = -((k(ind+1)+k(ind))./2).*...
    t.dt/N.dz^2./rho(ind)./cp(ind);

% Boundary conditions =================================================== %
switch lower(T.ubound)
    case {'direchlet','const'}
        diag(1,2)       =   1;
    case {'neumann','flux'}
        diag(1,2)    =   1 + ...
            ( (k(1)+k(2))./2 + (k(1)+k(2))./2 ).*...
            t.dt/N.dz^2./rho(1)./cp(1);
        diag(2,3)  =   ...
            -( (k(1)+k(2))./2 + (k(1)+k(2))./2 ).*...
            t.dt/N.dz^2./rho(1)./cp(1);        
end
switch lower(T.lbound)
    case {'direchlet','const'}
        diag(N.nz,2)    =   1;
    case {'neumann','flux'}
        diag(N.nz-1,1)  =   ...
            -( (k(N.nz-1)+k(N.nz))./2 + (k(N.nz-1)+k(N.nz))./2 ).*...
            t.dt/N.dz^2./rho(N.nz)./cp(N.nz);
        diag(N.nz,2)    =   1 + ...
            ( (k(N.nz-1)+k(N.nz))./2 + (k(N.nz-1)+k(N.nz))./2 ).*...
            t.dt/N.dz^2./rho(N.nz)./cp(N.nz);
    otherwise
        error('Boundary condition not defined!')
end

A       = spdiags(diag,[-1,0,1],N.nz,N.nz);

% Aenderung der rechten Seite durch andere Randbedingungen -------------- %
% switch lower(T.ubound)
%     case {'neumann','flux'}
%         T0(1)   = T0(1) + 2*s*N.dz*c1;
% end
switch lower(T.ubound)
    case {'neumann','flux'}
        T0(1)    =   T0(1) - ...
            2*T.utbf*( (k(2) + k(1))./2 ).*...
            t.dt/N.dz./rho(1)./cp(1);
end
switch lower(T.lbound)
    case {'neumann','flux'}
        T0(N.nz)    =   T0(N.nz) + ...
            2*T.ltbf*( (k(N.nz-1) + k(N.nz))./2 ).*...
            t.dt/N.dz./rho(N.nz)./cp(N.nz);
end

rhs     =   T0 + H.*t.dt./cp;

% Compute new temperature
T.T     =   A\rhs;

end
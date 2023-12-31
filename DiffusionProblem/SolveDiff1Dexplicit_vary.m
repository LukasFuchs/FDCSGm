function [T] = SolveDiff1Dexplicit_vary(N,T,Py,t)

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

ind     =   2:N.nz-1;

% Boundary conditions =================================================== %
switch lower(T.ubound)
    case {'direchlet','const'}
        T.T(1)      =   T0(1);         
    case {'neumann','flux'}
        kB          =   (k(2)+k(1))/2;
        kA          =   (k(1)+k(1))/2;
        
        a           =   (t.dt*(kA+kB))./(N.dz^2.*rho(1).*cp(1));
        b           =   1 - (t.dt.*(kB + kA))./(N.dz^2.*rho(1).*cp(1));
        c           =   -(kA*t.dt*2*T.utbf)/(N.dz*rho(1)*cp(1));
        
        T.T(1)      =   a*T0(2) + b*T0(1) + c + ...
                            H(1).*t.dt./cp(1);
end
switch lower(T.lbound)
    case {'direchlet','const'}
        T.T(N.nz)   =   T0(N.nz);
    case {'neumann','flux'}
        kB          =   (k(N.nz)+k(N.nz))/2;
        kA          =   (k(N.nz)+k(N.nz-1))/2;
                
        a           =   (t.dt*(kA+kB))./(N.dz^2.*rho(N.nz).*cp(N.nz));
        b           =   1 - (t.dt.*(kA + kB))./(N.dz^2.*rho(N.nz).*cp(N.nz));
        c           =   (kB*t.dt*2*T.ltbf)/(N.dz*rho(N.nz)*cp(N.nz));
        
        T.T(N.nz)   =   a*T0(N.nz-1) + b*T0(N.nz) + c + ...
                            H(N.nz).*t.dt./cp(N.nz);
    otherwise
        error('Boundary condition not defined!')
end

kA      =   (k(ind-1) + k(ind))/2;
kB      =   (k(ind) + k(ind+1))/2;

a       =   (kB.*t.dt)./(N.dz^2.*rho(ind).*cp(ind));

b       =   1 - (t.dt.*(kA + kB))./(N.dz^2.*rho(ind).*cp(ind));

c       =   (kA.*t.dt)./(N.dz^2.*rho(ind).*cp(ind));

T.T(ind)    =   a.*T0(ind+1) + b.*T0(ind) + c.*T0(ind-1) + ... 
                    H(ind).*t.dt./cp(ind);

end
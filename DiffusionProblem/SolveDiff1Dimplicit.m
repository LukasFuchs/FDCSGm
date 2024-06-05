function [T1,beenhere,A] = ...
    SolveHeat1Dimplicit(T0,nx,kappa,dt,dx,bound,A,beenhere,Q,rho,cp)
% Erstllung der Koeffizientenmatrix fuer implizites Problem ------------- %

if nargin < 9
    % Keine Waermequelle gegeben
    Q   = zeros(nx,1);      % [Q] = W/m^3; [Q] = [rho*H]
    rho = 1;
    cp  = 1;
end

s           = kappa*dt/dx^2;

switch bound
    % Definier Fluss an den Raendern
    case 'Neumann'
        c1  = 0;    c2  = 0;
end

% if beenhere == 0
ind         = 2:nx-1;

diag        = zeros(nx,3);

diag(ind-1,1)   = -s;
diag(ind,2)     = (1+2*s);
diag(ind+1,3)   = -s;

switch lower(bound)
    case {'direchlet','const'}
        diag(1,2)   = 1;
        diag(nx,2)  = 1;
    case {'neumann','flux'}
        % Linke Seite
        diag(1,2)   = (1+2*s);
        diag(2,3)   = -2*s;
        % Rechte Seite
        diag(nx-1,1)= -2*s;
        diag(nx,2)  = (1+2*s);
    otherwise
        error('Boundary condition not defined!')
end

A       = spdiags(diag,[-1,0,1],nx,nx);
% spy(A);
%     beenhere    = beenhere + 1;
% end

% Aenderung der rechten Seite durch andere Randbedingungen -------------- %
switch lower(bound)
    case {'neumann','flux'}
        T0(1)   = T0(1) + 2*s*dx*c1;
        T0(nx)  = T0(nx) + 2*s*dx*c2;
end

rhs     = T0 + Q.*dt/rho/cp;

% Compute new temperature
T1      = A\rhs;

end
function [F] = SolvePoisson2Dvaryk(Q,N,B,k)
% Function to solve 2D heat diffusion equation using the explicit finite
% difference scheme
% [Q] = W/m^3
% ----------------------------------------------------------------------- %

%% Define Coeffizients for Matrix A and boundary conditions ------------- %

% Erstellung der durchlaufenden Indizes ----------------------------- %
Number  = zeros(N.nz,N.nx);
num     = 1;
for i=1:N.nz
    for j=1:N.nx
        Number(i,j) = num;
        num = num+1;
    end
end

% Setup coefficient matrix A ---------------------------------------- %
a       = 1/N.dz^2;
b       = 1/N.dx^2;
% Reshape k into a vector so that the index fits later in the coefficient
% calculations
k       = reshape(k',[N.nx*N.nz,1]);

% Define diagonals for matrix
diag    = zeros(N.nx*N.nz,5);

% Inner index ------------------------------------------------------- %
iind   = Number(2:(N.nz-1),2:(N.nx-1));

diag(iind-N.nx,1)   = a.*(k(iind-N.nx,1)+k(iind,1))./2;
diag(iind-1,2)      = b.*(k(iind,1)+k(iind-1,1))./2;
diag(iind,3)        = ...
    -b.*((k(iind+1,1)+k(iind,1))./2+(k(iind,1)+k(iind-1,1))./2) ...
    - a.*((k(iind+N.nx,1)+k(iind,1))./2+(k(iind-N.nx,1)+k(iind-1,1))./2);
diag(iind+1,4)      = b.*(k(iind+1,1)+k(iind,1))./2;
diag(iind+N.nx,5)   = a.*(k(iind+N.nx,1)+k(iind,1))./2;

% Top and bottom index and boundary conditions ---------------------- %
switch B.ttbc
    case {'const','dirichlet'}
        tind                = Number(1,:);
        diag(tind,3)        = 1;
    case {'flux','neumann'}
        tind                = Number(1,1:N.nx);
        diag(tind,3)        = -1;
        diag(tind+N.nx,1)   = 1;
    otherwise
        error('Thermal boundary condition not defined!')
end
switch B.btbc
    case {'const','dirichlet'}
        bind                = Number(N.nz,:);
        diag(bind,3)        = 1;
    case {'flux','neumann'}
        bind                = Number(N.nz,1:N.nx);
        diag(bind-N.nx,1)   = 1;
        diag(bind,3)        = -1;
    otherwise
        error('Thermal boundary condition not defined!')
end

% Left and right index and boundary conditions ---------------------- %
lind    = Number(2:(N.nz1),1);
rind    = Number(2:(N.nz1),N.nx);

switch B.ltbc
    case {'const','dirichlet'}
        diag(lind,3)    = 1;
    case {'flux','neumann'}
        diag(lind,3)        = -1;
        diag(lind+1,4)      = 1;
    otherwise
        error('Thermal boundary condition not defined')
end
switch B.rtbc
    case {'const','dirichlet'}
        diag(rind,3)    = 1;
    case {'flux','neumann'}
        diag(rind-1,2)  = 1;
        diag(rind,3)    = -1;
    otherwise
        error('Thermal boundary condition not defined')
end

A       = spdiags(diag,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);

%% Solve system of equations -------------------------------------------- %
rhs    = -1.*reshape(Q',[N.nx*N.nz,1]);

switch B.btbc
    case {'const','dirichlet'}
        rhs(bind)   =   B.bhf;
    case {'flux','neumann'}
        %i=1
        rhs(bind)   = B.bhf*N.dz./...
            ((k(bind)+k(bind-N.nx))./2);
end
switch B.ttbc
    case {'const','dirichlet'}
        rhs(tind)   =   B.thf;
    case {'flux','neumann'}
        %i=nz
        rhs(tind)   = B.thf*N.dz./...
            ((k(tind+N.nx)+k(tind))./2);
end
switch B.ltbc
    case {'flux','neumann'}
        rhs(lind)   = B.lhf*N.dx./...
            ((k(lind)+k(lind+1))./2);
end
switch B.rtbc
    case {'flux','neumann'}
        rhs(rind)   = B.rhf*N.dx./...
            ((k(rind)+k(rind+1))./2);
end

F      = A\rhs;

F      = reshape(F,[N.nx,N.nz])';

end
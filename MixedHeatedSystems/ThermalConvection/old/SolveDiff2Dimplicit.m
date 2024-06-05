function [T1,N] = SolveDiff2Dimplicit(T0,Q,rho,dt,Py,N,B)
% Function to solve 2D heat diffusion equation using the explicit finite
% difference scheme
% assuming constant k, rho, cp
% dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
% ----------------------------------------------------------------------- %


%% Define Coeffizients for Matrix A and boundary conditions ------------- %
% Setup coefficient matrix A -------------------------------------------- %
sx      =   Py.kappa*dt/N.dx^2;
sz      =   Py.kappa*dt/N.dz^2;

% if N.beenhere == 0
% Erstellung der durchlaufenden Indizes ----------------------------- %
Number  =   zeros(N.nz,N.nx);
num     =   1;
for i=1:N.nz
    for j=1:N.nx
        Number(i,j) =   num;
        num = num+1;
    end
end

% Define diagonals for matrix
diag    =   zeros(N.nx*N.nz,5);

% Inner index ------------------------------------------------------- %
iind    =    Number(2:(N.nz-1),2:(N.nx-1));

diag(iind-N.nx,1)   =   -sz;
diag(iind-1,2)      =   -sx;
diag(iind,3)        =   1+2*sx+2*sz;
diag(iind+1,4)      =   -sx;
diag(iind+N.nx,5)   =   -sz;

% Top and bottom index and boundary conditions ---------------------- %
switch lower(B.btbc)
    case {'const','dirichlet'}
        bind                =   Number(N.nz,:);
        diag(bind,3)        =   1;
    case {'flux','neumann'}
        bind                =   Number(N.nz,2:(N.nx-1));
        diag(bind-N.nx,1)   =   -2*sz;
        diag(bind-1,2)      =   -sx;
        diag(bind,3)        =   1+2*sx+2*sz;
        diag(bind+1,4)      =   -sx;
        
        diag(N.nx*(N.nz-1)+1,3) =   1; % Left Corner grid point
        diag(N.nx*N.nz,3)       =   1; % Right Corner grid point
    otherwise
        error('Thermal boundary condition not defined!')
end
switch lower(B.ttbc)
    case {'const','dirichlet'}
        tind                =   Number(1,:);
        diag(tind,3)        =   1;
    case {'flux','neumann'}
        tind                =   Number(1,2:(N.nx-1));
        diag(tind-1,2)      =   -sx;
        diag(tind,3)        =   1+2*sx+2*sz;
        diag(tind+1,4)      =   -sx;
        diag(tind+N.nx,5)   =   -2*sz;
        
        diag(1,3)           =   1; % Left Corner grid point
        diag(N.nx,3)        =   1; % Right Corner grid point
    otherwise
        error('Thermal boundary condition not defined!')
end

% Left and right index and boundary conditions ---------------------- %
lind    = Number(2:(N.nz-1),1);
rind    = Number(2:(N.nz-1),N.nx);

switch lower(B.ltbc)
    case {'const','dirichlet'}
        diag(lind,3)        = 1;
    case {'flux','neumann'}
        diag(lind-N.nx,1)   = -sz;
        diag(lind,3)        = 1+2*sx+2*sz;
        diag(lind+1,4)      = -2*sx;
        diag(lind+N.nx,5)   = -sz;
    otherwise
        error('Thermal boundary condition not defined')
end
switch lower(B.rtbc)
    case {'const','dirichlet'}
        diag(rind,3)        = 1;
    case {'flux','neumann'}
        diag(rind-N.nx,1)   = -sz;
        diag(rind-1,2)      = -2*sx;
        diag(rind,3)        = 1+2*sx+2*sz;
        diag(rind+N.nx,5)   = -sz;
    otherwise
        error('Thermal boundary condition not defined')
end

N.A         = spdiags(diag,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);

%     N.beenhere  = N.beenhere + 1;
% end

%% Solve system of equations -------------------------------------------- %
switch lower(B.btbc)
    case {'flux','neumann'}
        T0(N.nz,2:(N.nx-1)) = T0(N.nz,2:(N.nx-1)) - 2*sz*N.dz*B.bhf;
    case {'const','dirichlet'}
        Q(N.nz,:)           = 0;
end
switch lower(B.ttbc)
    case {'flux','neumann'}
        T0(1,2:(N.nx-1))    = T0(1,2:(N.nx-1)) - 2*sz*N.dz*B.thf;
    case {'const','dirichlet'}
        Q(1,:)              = 0;
end
switch lower(B.ltbc)
    case {'flux','neumann'}
        T0(2:(N.nz-1),1)    = T0(2:(N.nz-1),1) + 2*sx*N.dx*B.lhf;
    case {'const','dirichlet'}
        Q(2:N.nz-1,1)       = 0;
end
switch lower(B.rtbc)
    case {'flux','neumann'}
        T0(2:(N.nz-1),N.nx) = T0(2:(N.nz-1),N.nx) + 2*sx*N.dx*B.rhf;
    case {'const','dirichlet'}
        Q(2:N.nz-1,N.nx)    = 0;
end

rhsT    = reshape(T0',[N.nx*N.nz,1]);
rhsQ    = reshape(Q',[N.nx*N.nz,1]).*dt./...
            reshape(rho',[N.nx*N.nz,1])./Py.cp;

rhs     = rhsT + rhsQ;

T1      = N.A\rhs;

T1      = reshape(T1,[N.nx,N.nz])';

switch lower(B.ttbc)
    case {'flux','neumann'}
        T1(1,1)         = T1(1,2);
        T1(1,N.nx)      = T1(1,N.nx-1);
end
switch lower(B.btbc)
    case {'flux','neumann'}
        T1(N.nz,1)       = T1(N.nz,2);
        T1(N.nz,N.nx)    = T1(N.nz,N.nx-1);
end

end
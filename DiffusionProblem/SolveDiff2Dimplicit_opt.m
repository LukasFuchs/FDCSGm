function To = SolveDiff2Dimplicit_opt(T0,Q,rho,dt,P,N,B)
% Function to solve 2D heat diffusion equation using the explicit finite
% difference scheme
% assuming constant k, rho, cp
% dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
% ----------------------------------------------------------------------- %

%% Define Coeffizients for Matrix A and boundary conditions ------------- %
To      = zeros(N.nz,N.nx);

% Setup coefficient matrix A -------------------------------------------- %
sx      =   P.kappa*dt/N.dx^2;
sz      =   P.kappa*dt/N.dz^2;


% Erstellung der durchlaufenden Indizes --------------------------------- %
neqx    =   N.nx;
neqz    =   N.nz;
switch B.ltbc
    case {'const','dirichlet'}
        neqx    =   neqx - 1;
end
switch B.rtbc
    case {'const','dirichlet'}
        neqx    =   neqx - 1;
end
switch B.ttbc
    case {'const','dirichlet'}
        neqz    =   neqz - 1;
end
switch B.btbc
    case {'const','dirichlet'}
        neqz    =   neqz - 1;
end

neq         =   neqx*neqz;

% if N.beenhere == 0
N.Number    =   zeros(N.nz,N.nx);
N.Number2   =   zeros(N.nz,N.nx);

num     =   1;
for j=2:(N.nx-1)
    for i=2:(N.nz-1)
        N.Number2(i,j)    =   num;
        num = num+1;
    end
end
num     =   1;
for j=1:N.nx
    for i=1:N.nz
        N.Number(i,j)     =   num;
        num = num+1;
    end
end

% Index der inneren Gitterpunkte ------------------------------------ %
N.iind              =   N.Number(2:(N.nz-1),2:(N.nx-1));                            %
%   Example index for a 6 x 6 FD-Problem
% Boundary indices for the global grid ------------------------------ %             %   iind
N.tind              =   N.Number(N.nz,3:N.nx-2);                                    %   ====================================
N.tlind             =   N.Number(N.nz,2);                                           %               j ->            Bottom
N.trind             =   N.Number(N.nz,N.nx-1);                                      %       1   7   13  19  25  31
%    i  2   8   14  20  26  32
N.lind              =   N.Number(3:(N.nz-2),1);                                     %    |  3   9   15  21  27  33
N.lbind             =   N.Number(2,1);                                              %    v  4  10   16  22  28  34
N.ltind             =   N.Number(N.nz-1,1);                                         %       5  11   17  23  29  35
%       6  12   18  24  30  36
N.rind              =   N.Number(3:(N.nz-2),N.nx);                                  %                                 Top
N.rbind             =   N.Number(2,N.nx);                                           %
N.rtind             =   N.Number(N.nz-1,N.nx);                                      %
%
N.bind              =   N.Number(1,3:N.nx-2);                                       %
N.blind             =   N.Number(1,2);                                              %
N.brind             =   N.Number(1,N.nx-1);                                         %

% Boundary indices for the inner grid in the matrix ----------------- %             %       iind2
N.tind2             =   N.Number2(N.nz-1,3:N.nx-2);                                 %       =====================
N.tlind2            =   N.Number2(N.nz-1,2);                                        %               j ->
N.trind2            =   N.Number2(N.nz-1,N.nx-1);                                   %       0   0   0   0   0   0
%    i  0   1   5   9   13  0
N.lind2             =   N.Number2(3:N.nz-2,2);                                      %    |  0   2   6   10  14  0
N.rind2             =   N.Number2(3:N.nz-2,N.nx-1);                                 %    v  0   3   7   11  15  0
%       0   4   8   12  16  0
N.bind2             =   N.Number2(2,3:N.nx-2);                                      %       0   0   0   0   0   0
N.blind2            =   N.Number2(2,2);                                             %
N.brind2            =   N.Number2(2,N.nx-1);                                        %
%

% Define diagonals for matrix --------------------------------------- %             %    Coefficient matrix - for Direchlet boundary conditions
diag2               =   zeros(neq,5);                                               %   ===================================================================
%   * -> zeros for boundary coefficient notes
diag2(1:(end-neqz),1)       =   -sx;            % W                                 %   ( C   N   ... ... E   ... ... ... ... ... ... ... ... ... ... ... )
diag2(1:(end-1),2)          =   -sz;            % S                                 %   ( S   C   N   ... ... E   ... ... ... ... ... ... ... ... ... ... )-
diag2(:,3)                  =   1+2*sx+2*sz;    % C                                 %   ( ... S   C   N   ... ... E   ... ... ... ... ... ... ... ... ... )
diag2(2:end,4)              =   -sz;            % N                                 %   ( ... ... S   C   *   ... ... E   ... ... ... ... ... ... ... ... )
diag2((neqz+1):end,5)       =   -sx;            % E                                 %   ( W   ... ... *   C   N   ... ... E   ... ... ... ... ... ... ... )
%   ( ... W   ... ... S   C   N   ... ... E   ... ... ... ... ... ... )
% Remove coefficients at the boundaries for diagonals S and N                       %   ( ... ... W   ... ... S   C   N   ... ... E   ... ... ... ... ... )
diag2(N.bind2-1,2)          = 0;                                                    %   ( ... ... ... W   ... ... S   C   *   ... ... E   ... ... ... ... )
diag2(N.brind2-1,2)         = 0;                                                    %   ( ... ... ... ... W   ... ... *   C   N   ... ... E   ... ... ... )
diag2(N.tind2+1,4)          = 0;                                                    %   ( ... ... ... ... ... W   ... ... S   C   N   ... ... E   ... ... )
diag2(N.tlind2+1,4)         = 0;                                                    %   ( ... ... ... ... ... ... W   ... ... S   C   N   ... ... E   ... )
%   ( ... ... ... ... ... ... ... W   ... ... S   C   *   ... ... E   )
N.A     = spdiags(diag2,[-neqz,-1,0,1,neqz],...                                     %   ( ... ... ... ... ... ... ... ... W   ... ... *   C   N   ... ... )
    neq,neq);                                                                       %   ( ... ... ... ... ... ... ... ... ... W   ... ... S   C   N   ... )
%   ( ... ... ... ... ... ... ... ... ... ... W   ... ... S   C   N   )
%     N.beenhere  = N.beenhere + 1;                                                       %   ( ... ... ... ... ... ... ... ... ... ... ... W   ... ... S   C   )
%
% end                                                                                     %

rhsT    = reshape(T0(N.iind),[neq,1]);                                                  %   RHS - for Direchlet boundary conditions
rhsQ    = reshape(Q(N.iind),[neq,1]).*dt./...
    reshape(rho(N.iind),[neq,1])./P.cp;                                         %   ======================================================
%
% Left and right index and boundary conditions -------------------------- %             %   ( T0(iind2(1)  + sz*T0(iind1(7)) + sx*T0(iind1(2))   )
switch B.ltbc                                                                           %   ( T0(iind2(2)  + sx*T0(iind1(3))                     )
    case {'const','dirichlet'}                                                          %   ( T0(iind2(3)  + sx*T0(iind1(4))                     )
        % j = 1                                                                         %   ( T0(iind2(4)  + sz*T0(iind1(12)) + sx*T0(iind1(5))  )
        % T(i,j) + sx*T(i,j-1) + sz*T(i-1,j)/sz*T(i+1,j)                                %   ( T0(iind2(5)  + sz*T0(iind1(13))                    )
        rhsT(N.blind2)     =   rhsT(N.blind2) + sz*T0(N.blind) + sx*T0(N.lbind);        %   ( T0(iind2(6)                                        )
        rhsT(N.lind2)      =   rhsT(N.lind2) +  sx*T0(N.lind);                          %   ( T0(iind2(7)                                        )
        rhsT(N.tlind2)     =   rhsT(N.tlind2) + sz*T0(N.tlind) + sx*T0(N.ltind);        %   ( T0(iind2(8)  + sz*T0(iind1(18))                    )
    otherwise                                                                           %   ( T0(iind2(9)  + sz*T0(iind1(19))                    )
        error('Thermal boundary condition not defined')                                 %   ( T0(iind2(10)                                       )
end                                                                                     %   ( T0(iind2(11)                                       )
switch B.btbc                                                                           %   ( T0(iind2(12) + sz*T0(iind1(24))                    )
    case {'const','dirichlet'}                                                          %   ( T0(iind2(13) + sz*T0(iind1(25)) + sx*T0(iind1(32)) )
        % i = 1                                                                         %   ( T0(iind2(14) + sx*T0(iind1(33))                    )
        % T(i,j) + sz*T(i-1,j) + sx*T(i,j-1)/sx*T(i,j+1)                                %   ( T0(iind2(15) + sx*T0(iind1(34))                    )
        %   ( T0(iind2(16) + sz*T0(iind1(30)) + sx*T0(iind1(35)) )
        rhsT(N.bind2)       =   rhsT(N.bind2) + sz*T0(N.bind)';
        
    otherwise
        error('Thermal boundary condition not defined!')
end
switch B.ttbc
    case {'const','dirichlet'}
        % i = nz
        % T(i,j) + sz*T(i+1,j) + sx*T(i,j-1)/sx*T(i,j+1)
        rhsT(N.tind2)       =   rhsT(N.tind2) + sz*T0(N.tind)';
        
    otherwise
        error('Thermal boundary condition not defined!')
end
switch B.rtbc
    case {'const','dirichlet'}
        % j = nx
        % T(i,j) + sx*T(i,j+1) + sz*T(i-1,j)/sz*T(i+1,j)
        rhsT(N.brind2)      =   rhsT(N.brind2) + sz*T0(N.brind) + sx*T0(N.rbind);
        rhsT(N.rind2)       =   rhsT(N.rind2) +  sx*T0(N.rind);
        rhsT(N.trind2)      =   rhsT(N.trind2) + sz*T0(N.trind) + sx*T0(N.rtind);
    otherwise
        error('Thermal boundary condition not defined')
end

%% Solve system of equations -------------------------------------------- %
% switch B.ttbc
% %     case {'flux','neumann'}
% %         T0(N.nz,2:(N.nx-1)) = T0(N.nz,2:(N.nx-1)) - 2*sz*N.dz*B.thf;
%     case {'const','dirichlet'}
%         Q(N.nz,:)           = 0;
% end
% switch B.btbc
% %     case {'flux','neumann'}
% %         T0(1,2:(N.nx-1))    = T0(1,2:(N.nx-1)) - 2*sz*N.dz*B.bhf;
%     case {'const','dirichlet'}
%         Q(1,:)              = 0;
% end
% switch B.ltbc
% %     case {'flux','neumann'}
% %         T0(2:(N.nz-1),1)    = T0(2:(N.nz-1),1) + 2*sx*N.dx*B.lhf;
%     case {'const','dirichlet'}
%         Q(2:N.nz-1,1)       = 0;
% end
% switch B.rtbc
% %     case {'flux','neumann'}
% %         T0(2:(N.nz-1),N.nx) = T0(2:(N.nz-1),N.nx) + 2*sx*N.dx*B.rhf;
%     case {'const','dirichlet'}
%         Q(2:N.nz-1,N.nx)    = 0;
% end

% rhsT    = reshape(T0',[N.nx*N.nz,1]);
% rhsQ    = reshape(Q',[N.nx*N.nz,1]).*dt/P.rho/P.cp;

rhs         = rhsT + rhsQ;

% tic;
T1          = N.A\rhs;
% tend        = toc;
% tend

To(N.iind)        = reshape(T1,[neqz,neqx]);

To(~N.Number2)    = T0(~N.Number2);
% keyboard

% switch lower(B.btbc)
%     case {'flux','neumann'}
%         T1(1,1)         = T1(1,2);
%         T1(1,N.nx)      = T1(1,N.nx-1);
% end
% switch lower(B.ttbc)
%     case {'flux','neumann'}
%         T1(N.nz,1)       = T1(N.nz,2);
%         T1(N.nz,N.nx)    = T1(N.nz,N.nx-1);
% end

end
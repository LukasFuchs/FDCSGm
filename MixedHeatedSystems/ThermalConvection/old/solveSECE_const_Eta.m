function [D,A] = solveSECE_const_Eta(D,Py,N,B,A)
% Funkion
% ----------------------------------------------------------------------- %
% Scaling parameters for linear equations solver ------------------------ %
N.Kbond     =   4*Py.eta0/(N.dx-N.dz)^2;
N.Kcont     =   2*Py.eta0/(N.dx-N.dz);

if ~(B.beenhere == 1)
    % Rescaling of boundary conditions ---------------------------------- %
    D.vxi       =   D.vxi.*N.Kbond;
    D.vzi       =   D.vzi.*N.Kbond;
    D.Pi(B.B0c) =   D.Pi(B.B0c)./N.Kcont.*N.Kbond;
end

% SYSTEM OF EQUATIONS --------------------------------------------------- %
N.np        =   3;                  % NUMBER OF PARAMETERS AT A GRID POINT
N.NN        =   N.np*N.nx*N.nz;     % TOTAL NUMBER OF EQUATIONS

%%   BUILD UP THE RHS --------------------------------------------------- %
rhs     = zeros(N.NN,1);

for j = 1:N.nx
    for i = 1:N.nz
        G = N.np*(N.nz*(j-1)+i);
        rhs(G-2,1)  =   D.Pi(i,j);
        rhs(G-1,1)  =   D.vxi(i,j);
        rhs(G,1)    =   D.vzi(i,j);
        if (j>1)&&(i>1)&&(j<(N.nx1))&&(i<N.nz)
            rhs(G,1)    =   D.vzi(i,j) + ...
                Py.g*(D.rho(i,j) + D.rho(i,j+1)) / 2;
        end
    end
end

%%  BUILDING COEFFICIENT MATRIX ----------------------------------------- %
if ~(B.beenhere == 1)
    A       =   sparbuistag_constEta(N,Py.eta0,B);
end

%% SOLVING SYSTEM OF EQUATIONS ------------------------------------------ %
b       =   A\rhs;

%% REALLOCATE THE VALUES TO THE GRID ------------------------------------ %
for j = 1:N.nx
    for i = 1:N.nz
        G   =   N.np*(N.nz*(j-1)+i);
        
        D.P(i,j)    =     b(G-2,1)*N.Kcont;
        D.vx(i,j)   =     b(G-1,1);
        D.vz(i,j)   =     b(G,1);
    end
end

end

function [A] = sparbuistag_constEta(N,eta,B)

%% SETTING INDEXES FOR THE DIAGONAL NUMBERS ----------------------------- %
sgp         =   7;              % NUMBER OF STENCIL GRID POINTS
np          =   3;              % NUMBER OF PARAMETERS PER GRID POINT
NN          =   3*N.nx*N.nz;    % TOTAL NUMBER OF EQUATIONS
NNZD        =   3*sgp+3;        % NUMBER OF NON-ZERO DIAGONALS

% ---------- Setting diagonal vector and coefficient bandmatrix --------- %
di          = zeros(NNZD,1);
ABC         = zeros(NN,NNZD);

%% DIAGONAL INDEX ------------------------------------------------------- %
[di]        = diagnumber_constEta(di,np,N.nz);

%% BUILDING THE SPARSE MATRIX ------------------------------------------- %
% COLLECTING COEFFICIENTS FROM THE GRID --------------------------------- %
for j = 1:N.nx
    for i = 1:N.nz
        G = np*( N.nz*(j-1) + i);
        
        %% GHOST NODE COEFFICIENTS -------------------------------------- %
        % ----------------------- X-STOKES EQUATION --------------------- %
        if (i==N.nz)
            ABC(G-1+di(10),10)  =   N.Kbond;
        end
        % ----------------------- Z-STOKES EQUATION --------------------- %
        if (j==N.nx)
            ABC(G+di(10),10)    =   N.Kbond;
        end
        % ---------------------- CONTINUUM EQUATION --------------------- %
        if(j==1)||(i==1)
            ABC(G-2+di(10),10)  =   N.Kbond;
        end
        
        %% X-STOKES EQUATION -------------------------------------------- %
        ABC     =   xstokes_const_Eta(ABC,G,di,i,j,N,eta,B);
        %% Z-STOKES EQUATION -------------------------------------------- %
        ABC     =   zstokes_const_Eta(ABC,G,di,i,j,N,eta,B);
        %% CONTINUITY EQUATION ------------------------------------------ %
        ABC     =   contequat_const_Eta(ABC,G,di,i,j,N,B);
    end
end

A       = spdiags(ABC,di,NN,NN);
% figure(1)
% spy(A)
% keyboard
clear di ABC

end

function [ABC] = xstokes_const_Eta(ABC,G,di,i,j,N,eta,B)
% Function to collect the coefficients for the x-stokes equation using a
% staggered grid finite difference scheme with constant viscosity.
%
% ----------------------------------------------------------------------- %

%% BOUNDARY COEFFIECIENTS ----------------------------------------------- %
if (j==1)&&(i<N.nz)||(i==1)&&(j>1)&&(j<N.nx)||(j==N.nx)&&(i<N.nz)||...
        (i==(N.nz-1))&&(j>1)&&(j<N.nx)
    %% Left and right boundary ------------------------------------------ %
    % ------------- free slip or no slip -------------------------------- %
    % vx(i,j) = 0; or given velocity ------------------------------------ %
    ABC(G-1+di(10),10)      =   N.Kbond;                  % vx(i,j)
    
    %% UPPER BOUNDARY --------------------------------------------------- %
    switch B.tbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ----------- vx = 0: vx(i,j)-1/3*vx(i+1,j)=0 --------------- %
            if(i==1)&&(j>1)&&(j<N.nx)
                % South
                ABC(G-1+di(13),13)      = -N.Kbond/3;     % vx(j+1,i)
            end
        case 1
            % --------------------- FREE SLIP --------------------------- %
            % ----------- dvx/dz = 0; vx(i,j)-vx(i+1,j)=0 --------------- %
            if(i==1)&&(j>1)&&(j<N.nx)
                % South
                ABC(G-1+di(13),13)      = -N.Kbond;       % vx(j+1,i)
            end
    end
    %% LOWER BOUNDARY --------------------------------------------------- %
    switch B.bbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ----------- vx = 0: vx(i,j)-1/3*vx(i-1,j)=0 --------------- %
            if(i==(N.nz1))&&(j>1)&&(j<N.nx)
                % North
                ABC(G-1+di(7),7)      = -N.Kbond/3;     % vx(i-1,j)
            end
        case 1
            % --------------------- FREE SLIP --------------------------- %
            % ----------- dvx/dz = 0; vx(i,j)-vx(i-1,j)=0 --------------- %
            if(i==(N.nz1))&&(j>1)&&(j<N.nx)
                % North
                ABC(G-1+di(7),7)      = -N.Kbond;       % vx(i-1,j)
            end
    end
end

%% INNER GRID POINTS ---------------------------------------------------- %
if(i>1)&&(j>1)&&(j<N.nx)&&(i<(N.nz1))
    %% W
    ABC(G-1+di(4),4)    =   eta/N.dx^2;             % vx(i,j-1)
    %% N
    ABC(G-1+di(7),7)    =   eta/N.dz^2;             % vx(i-1,j)
    %% C
    ABC(G-1+di(10),10)  =   -2*eta*(1/N.dx^2 + 1/N.dz^2);   % vx(i,j)
    %% S
    ABC(G-1+di(12),12)  =   N.Kcont/N.dx;           % P(i+1,j)
    ABC(G-1+di(13),13)  =   eta/N.dz^2;             % vx(i+1,j)
    %% E
    ABC(G-1+di(19),19)  =   eta/N.dx^2;             % vx(i,j+1)
    %% SE
    ABC(G-1+di(21),21)  =   -N.Kcont/N.dx;          % P(i+1,j+1)
end
end

function [ABC] = zstokes_const_Eta(ABC,G,di,i,j,N,eta,B)
% Function to collect the coefficients for the z-stokes equation using a
% staggered grid finite difference scheme with constant viscosity.
%
% ----------------------------------------------------------------------- %

%% BOUNDARY COEFFIECIENTS ----------------------------------------------- %
if (j==1)||(i==1)&&(j>1)&&(j<(N.nx-1))||(j==(N.nx-1))||...
        (i==N.nz)&&(j>1)&&(j<(N.nx-1))
    %% Upper and lower boundary ----------------------------------------- %
    %  free, no slip or prescribed velocity ----------------------------- %
    %  In the case of fs/ns the velocity has to be zero! ---------------- %
    ABC(G+di(10),10)      = N.Kbond;              % vz(i,j)
    
    %% LEFT BOUNDARY CONDITION ------------------------------------------ %
    switch B.lbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ---------- vz = 0; vz(i,j) - vz(i,j+1)/3 = 0 -------------- %
            if (j==1)
                % East
                ABC(G+di(19),19)  =     -N.Kbond/3;     % vz(i,j+1)
            end
        case 1
            % --------------------- FREE SLIP --------------------------- %
            % --------- dvz/N.dx = 0; vz(i,j) - vz(i,j+1) = 0 ----------- %
            if (j==1)
                % East
                ABC(G+di(19),19)  =     -N.Kbond;       % vz(i,j+1)
            end
    end
    %% RIGHT BOUNDARY CONDITION ----------------------------------------- %
    switch B.rbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ------- vz = 0; vz(i,j) - 1/3*vz(i,j-1) = 0 --------------- %
            if (j==(N.nx-1))
                % West
                ABC(G+di(4),4)    =     -N.Kbond/3;     % vz(i,j-1)
            end
        case 1
            % -------------------- FREE SLIP ---------------------------- %
            % -------- dvz/dx = 0; vz(i,j) - vz(i,j-1) = 0 -------------- %
            if (j==(N.nx-1))
                ABC(G+di(4),4)    =     -N.Kbond;       % vz(j,i-1)
            end
    end
end

%% INNER GRID POINTS ---------------------------------------------------- %
if(j>1)&&(i>1)&&(j<(N.nx-1))&&(i<N.nz)
    %% W
    ABC(G+di(4),4)      =   eta/N.dx^2;                 % vz(i,j-1)
    %% N
    ABC(G+di(7),7)      =   eta/N.dz^2;                 % Vz(i-1,j)
    %% C
    ABC(G+di(10),10)    =   -2*eta*(1/N.dx^2 + 1/N.dz^2);    % Vz(i,j)
    %% S
    ABC(G+di(13),13)    =   eta/N.dz^2;                 % Vz(i+1,j)
    %% E
    ABC(G+di(17),17)    =   N.Kcont/N.dz;               % P(i,j+1)
    ABC(G+di(19),19)    =   eta/N.dx^2;                 % vz(i,j+1)
    %% SE
    ABC(G+di(20),20)    =   -N.Kcont/N.dz;              % P(i+1,j+1)
end
end

function [ABC] = contequat_const_Eta(ABC,G,di,i,j,N,B)


%% BOUNDARY COEFFIECIENTS ----------------------------------------------- %
if(i>1)&&(j>1)
    switch B.rbc
        case {0,1}
            if (j==2)&&(i==2) || (j==2)&&(i==N.nz) || (j==3)&&(i==2) || ...
                    (j==N.nx)&&(i==2) || (j==N.nx)&&(i==N.nz)
                %% ONE INTERIOR CELL ------------------------------------ %
                ABC(G-2+di(10),10)    =     N.Kbond;        % P(i,j)
                %% UPPER AND LOWER LEFT CORNER -------------------------- %
                % --------------- dP/dx=0 => P(i,j)-P(i,j+1) = 0 -------- %
                if (i==2)&&(j==2) || (i==N.nz)&&(j==2)
                    % East
                    ABC(G-2+di(19),19)=     -N.Kbond;       % P(i,j+1)
                end
                %% UPPER AND LOWER RIGHT CORNER ------------------------- %
                % --------------- dP/dx=0 => P(j,i-1)-P(j,i) = 0 -------- %
                if (i==2)&&(j==N.nx) || (i==N.nz)&&(j==N.nx)
                    % West
                    ABC(G-2+di(4),4)  =     -N.Kbond;       % P(i,j-1)
                end
            else
                %% INNER GRID POINTS ------------------------------------- %
                %% NW
                ABC(G-2+di(2),2)    =   -N.Kcont/N.dx;      % vx(i-1,j-1)
                ABC(G-2+di(3),3)    =   -N.Kcont/N.dz;      % vz(i-1,j-1)
                %% W
                ABC(G-2+di(6),6)    =   N.Kcont/N.dz;       % vz(i,j-1)
                %% N
                ABC(G-2+di(8),8)    =   N.Kcont/N.dx;       % vx(i-1,j)
            end
    end
end
end

function [di] = diagnumber_constEta(di,np,nz)
%% NW
di(1)       = -np*nz-np;
di(2)       = di(1)+1;
di(3)       = di(1)+2;
%% W
di(4)       = -np*nz;
di(5)       = di(4)+1;
di(6)       = di(4)+2;
%% N
di(7)       = -np;
di(8)       = di(7)+1;
di(9)       = di(7)+2;
%% C
di(10)      = 0;
di(11)      = di(10)+1;
di(12)      = di(10)+2;
%% S
di(13)      = np;
di(14)      = di(13)+1;
di(15)      = di(13)+2;
%% NE
di(16)      = np*nz-np;
di(17)      = di(16)+1;
di(18)      = di(16)+2;
%% E
di(19)      = np*nz;
di(20)      = di(19)+1;
di(21)      = di(19)+2;
%% SE
di(22)      = np*nz+np;
di(23)      = di(22)+1;
di(24)      = di(22)+2;
end



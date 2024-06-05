function [D,A] = solveSECE(D,Py,N,B)
% ----------------------------------------------------------------------- %
% Scaling parameters for linear equations solver ------------------------ %

N.Kbond     =   4*min(min(D.eta))/(N.dx-N.dz)^2;
N.Kcont     =   2*min(min(D.eta))/(N.dx-N.dz);

if ~(N.beenhere == 1)
    % Rescaling of boundary conditions ---------------------------------- %
    D.vxi       =   D.vxi.*N.Kbond;
    D.vzi       =   D.vzi.*N.Kbond;
    D.Pi(B.B0c) =   D.Pi(B.B0c)./N.Kcont.*N.Kbond;
end

% figure(1)
% subplot(1,3,1)
% plot(D.vzi(1:N.nz,1)./N.Kbond,'--')
% hold on
% % plot(D.vxi(1:N.nz,1)./N.Kbond,'-')
% set(gca,'yscale','log');
% subplot(2,3,2)
% % plot(D.vxi(1,1:N.nx)./N.Kbond,'-')
% hold on
% plot(D.vzi(1,1:N.nx)./N.Kbond,'--')
% set(gca,'yscale','log');
% subplot(2,3,5)
% % plot(D.vxi(N.nz1,1:N.nx)./N.Kbond,'-')
% hold on
% plot(D.vzi(N.nz1,1:N.nx)./N.Kbond,'--')
% set(gca,'yscale','log');
% subplot(1,3,3)
% plot(D.vzi(1:N.nz,N.nx1)./N.Kbond,'--')
% hold on
% % plot(D.vxi(1:N.nz,N.nx)./N.Kbond,'-')
% set(gca,'yscale','log');

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
A       =   sparbuistag(N,D.eta,B);

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
% figure(1)
% subplot(1,3,1)
% vxl = D.vx(1:N.nz1,1);
% vzl = D.vz(:,1) + (D.vz(:,1) - D.vz(:,2))/2;
% % plot(vxl,'.')
% hold on
% plot(vzl,'*')
% set(gca,'yscale','log');
% subplot(2,3,2)
% % vxtop = D.vx(1,:) + (D.vx(1,:) - D.vx(2,:))/2;
% vztop = D.vz(1,1:N.nx1);
% % plot(vxtop,'.')
% hold on
% plot(vztop,'*')
% set(gca,'yscale','log');
% subplot(2,3,5)
% % vxbot   = D.vx(N.nz1,:) + (D.vx(N.nz1,:) - D.vx(N.nz1-1,:))/2;
% vzbot   = D.vz(N.nz,1:N.nx1);
% % plot(vxbot,'.')
% hold on
% plot(vzbot,'*')
% set(gca,'yscale','log');
% subplot(1,3,3)
% % vxr = D.vx(1:N.nz1,N.nx);
% vzr = D.vz(:,N.nx1) + (D.vz(:,N.nx1) - D.vz(:,N.nx1-1))/2;
% % plot(vxr,'.')
% hold on
% plot(vzr,'*')
% set(gca,'yscale','log');

end

function [A,AA] = sparbuistag(N,eta,B)

%% SETTING INDEXES FOR THE DIAGONAL NUMBERS ----------------------------- %
sgp         =   3*4+2;          % NUMBER OF STENCIL GRID POINTS
np          =   3;              % NUMBER OF PARAMETERS PER GRID POINT
NN          =   3*N.nx*N.nz;    % TOTAL NUMBER OF EQUATIONS
NNZD        =   3*sgp;          % NUMBER OF NON-ZERO DIAGONALS

% ---------- Setting diagonal vector and coefficient bandmatrix --------- %
di          = zeros(NNZD,1);
ABC         = zeros(NN,NNZD);

%% DIAGONAL INDEX ------------------------------------------------------- %
[di]        = diagnumber(di,np,N.nz);

%% BUILDING THE SPARSE MATRIX ------------------------------------------- %
% COLLECTING COEFFICIENTS FROM THE GRID --------------------------------- %

for j = 1:N.nx
    for i = 1:N.nz
        G = np*( N.nz*(j-1) + i);
        
        %% GHOST NODE COEFFICIENTS -------------------------------------- %
        % ----------------------- X-STOKES EQUATION --------------------- %
        if (i==N.nz)
            ABC(G-1+di(19),19)  =   N.Kbond;
        end
        % ----------------------- Z-STOKES EQUATION --------------------- %
        if (j==N.nx)
            ABC(G+di(19),19)    =   N.Kbond;
        end
        % ---------------------- CONTINUUM EQUATION --------------------- %
        if(j==1)||(i==1)
            ABC(G-2+di(19),19)  =   N.Kbond;
        end
        
        %% X-STOKES EQUATION -------------------------------------------- %
        ABC     =   xstokes(ABC,G,di,i,j,N,eta,B);
        %% Z-STOKES EQUATION -------------------------------------------- %
        ABC     =   zstokes(ABC,G,di,i,j,N,eta,B);
        %% CONTINUITY EQUATION ------------------------------------------ %
        ABC     =   contequat(ABC,G,di,i,j,N,B);
    end
end

A       = spdiags(ABC,di,NN,NN);
% figure(1)
% spy(A)
% keyboard
clear di ABC

end

function [ABC] = xstokes(ABC,G,di,i,j,N,eta,B)
% Function to collect the coefficients for the x-stokes equation using a
% staggered grid finite difference scheme with variable viscosity.
%
% ----------------------------------------------------------------------- %

%% BOUNDARY COEFFIECIENTS ----------------------------------------------- %
% if (j==1)&&(i<N.nz)||(i==1)&&(j>1)&&(j<N.nx)||(j==N.nx)&&(i<N.nz)||...
%         (i==(N.nz-1))&&(j>1)&&(j<N.nx)
if (j==1)||(i==1)||(j==N.nx)||(i==(N.nz1))
    %% Left and right boundary ------------------------------------------ %
    % ------------- free, no slip and prescribed velocity --------------- %
    % vx(i,j) = 0; or given velocity ------------------------------------ %
    ABC(G-1+di(19),19)      =   N.Kbond;                  % vx(i,j)
    
    %% UPPER BOUNDARY --------------------------------------------------- %
    switch B.tbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ----------- vx = 0: 3*vx(i,j)/2-vx(i+1,j)/2=vxb ----------- %
            if(i==1)&&(j>1)&&(j<N.nx)
		ABC(G-1+di(19),19)      = 3*N.Kbond/2;      % vx(i,j)
                % South
                ABC(G-1+di(22),22)      = -N.Kbond/2;       % vx(i+1,j)
            end
        case 1
            % --------------------- FREE SLIP --------------------------- %
            % ----------- dvx/dz = 0; vx(i,j)-vx(i+1,j)=0 --------------- %
            if(i==1)&&(j>1)&&(j<N.nx)
                % South
                ABC(G-1+di(22),22)      = -N.Kbond;         % vx(i+1,j)
            end
    end
    %% LOWER BOUNDARY --------------------------------------------------- %
    switch B.bbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ----------- vx = 0: 3*vx(i,j)/2-vx(i-1,j)/2=vxb ----------- %
            if(i==(N.nz1))&&(j>1)&&(j<N.nx)
		ABC(G-1+di(19),19)      = 3*N.Kbond/2;      % vx(i,j)
                % North
                ABC(G-1+di(16),16)      = -N.Kbond/2;       % vx(i-1,j)
            end
        case 1
            % --------------------- FREE SLIP --------------------------- %
            % ----------- dvx/dz = 0; vx(i,j)-vx(i-1,j)=0 --------------- %
            if(i==(N.nz1))&&(j>1)&&(j<N.nx)
                % North
                ABC(G-1+di(16),16)      = -N.Kbond;       % vx(i-1,j)
            end
    end
end

%% INNER GRID POINTS ---------------------------------------------------- %
if(i>1)&&(j>1)&&(j<N.nx)&&(i<(N.nz1))
    C1      = -(eta(i,j-1)+eta(i+1,j-1)+2*(eta(i,j)+eta(i+1,j))+...
        eta(i,j+1)+eta(i+1,j+1))/2/N.dx^2;      % Mixed coefficients for vx(i,j)
    C2      = -(eta(i+1,j)+eta(i,j))/N.dz^2;    % Mixed coefficients for vx(i,j)
    
    %% W1
    ABC(G-1+di(4),4)      = 2/N.dx^2*(eta(i,j-1)+eta(i,j)+eta(i+1,j-1)+...
        eta(i+1,j))/4;                              % Vx
    ABC(G-1+di(5),5)      = eta(i,j)/N.dx/N.dz;     % Vz
    %% S1W
    ABC(G-1+di(8),8)      = -eta(i+1,j)/N.dx/N.dz;  % Vz
    %% N1
    ABC(G-1+di(16),16)     = eta(i,j)/N.dz^2;       % Vx
    %% Z
    ABC(G-1+di(19),19)     = C1+C2;                 % Vx
    ABC(G-1+di(20),20)     = -eta(i,j)/N.dx/N.dz;   % Vz
    %% S1
    ABC(G-1+di(21),21)     = N.Kcont/N.dx;          % Pressure at S1
    ABC(G-1+di(22),22)     = eta(i+1,j)/N.dz^2;     % Vx
    ABC(G-1+di(23),23)     = eta(i+1,j)/N.dx/N.dz;  % Vz
    %% E
    ABC(G-1+di(34),34)     = 1/N.dx^2*(eta(i,j)+eta(i,j+1)+eta(i+1,j)+...
        eta(i+1,j+1))/2;                            % Vx
    %% S1E1
    ABC(G-1+di(36),36)     = -N.Kcont/N.dx;         % Pressure at S1E1
end
end

function [ABC] = zstokes(ABC,G,di,i,j,N,eta,B)
% Function to collect the coefficients for the z-stokes equation using a
% staggered grid finite difference scheme with variable viscosity.
%
% ----------------------------------------------------------------------- %

%% BOUNDARY COEFFIECIENTS ----------------------------------------------- %
% if (j==1)||(i==1)&&(j>1)&&(j<(N.nx-1))||(j==(N.nx-1))||...
%         (i==N.nz)&&(j>1)&&(j<(N.nx-1))
if (j==1)||(i==1)||(j==(N.nx-1))||(i==N.nz)    
    %% Upper and lower boundary ----------------------------------------- %
    %  free, no slip or prescribed velocity ----------------------------- %
    %  In the case of fs/ns the velocity has to be zero! ---------------- %
    ABC(G+di(19),19)      = N.Kbond;              % vz(i,j)
    
    %% LEFT BOUNDARY CONDITION ------------------------------------------ %
    switch B.lbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ---------- vz = 0; 3*vz(i,j)/2 - vz(i,j+1)/2 = vzb -------- %
            if (j==1)&&(i>1)&&(i<(N.nz))
		ABC(G+di(19),19)  = 3*N.Kbond/2;    % vz(i,j)
                % East
                ABC(G+di(34),34)  = -N.Kbond/2;     % vz(i,j+1)
            end
        case 1
            % --------------------- FREE SLIP --------------------------- %
            % --------- dvz/N.dx = 0; vz(i,j) - vz(i,j+1) = 0 ----------- %
            if (j==1)&&(i>1)&&(i<(N.nz))
                % East
                ABC(G+di(34),34)  = -N.Kbond;     % vz(i,j+1)
            end
    end
    %% RIGHT BOUNDARY CONDITION ----------------------------------------- %
    switch B.rbc
        case 0
            % ---------------------- NO SLIP ---------------------------- %
            % ------- vz = 0; 3*vz(i,j)/2 - vz(i,j-1)/2 = 0 ------------- %
            if (j==(N.nx-1))&&(i>1)&&(i<(N.nz))
		ABC(G+di(19),19)  = 3*N.Kbond/2;    % vz(i,j)
                % West
                ABC(G+di(4),4)    = -N.Kbond/2;   % vz(i,j-1)
            end
        case 1
            % -------------------- FREE SLIP ---------------------------- %
            % -------- dvz/dx = 0; vz(i,j) - vz(i,j-1) = 0 -------------- %
            if (j==(N.nx-1))&&(i>1)&&(i<(N.nz))
                % West
                ABC(G+di(4),4)    = -N.Kbond;     % vz(j,i-1)
            end
    end
end

%% INNER GRID POINTS ---------------------------------------------------- %
if(j>1)&&(i>1)&&(j<(N.nx-1))&&(i<N.nz)
    C1      = -(eta(i-1,j)+eta(i-1,j+1)+2*(eta(i,j)+eta(i,j+1))+...
        eta(i+1,j)+eta(i+1,j+1))/2/N.dz^2;      % Mixed coefficient for vz(i,j)
    C2      = -(eta(i,j+1)+eta(i,j))/N.dx^2;    % Mixed coefficient for vz(i,j)
    
    %% W1
    ABC(G+di(4),4)       = eta(i,j)/N.dx^2;     % Vz
    %% N2
    ABC(G+di(15),15)     = eta(i,j)/N.dx/N.dz;  % Vx at N1
    %% N1
    ABC(G+di(16),16)     = (eta(i-1,j)+eta(i-1,j+1)+eta(i,j)+eta(i,j+1))/...
        2/N.dz^2;       % Vz
    ABC(G+di(18),18)     = -eta(i,j)/N.dx/N.dz; % Vx
    %% Z
    ABC(G+di(19),19)     = C1+C2;               % Vz
    %% S1
    ABC(G+di(22),22)     = (eta(i,j)+eta(i,j+1)+eta(i+1,j)+eta(i+1,j+1))/...
        2/N.dz^2;           % Vz
    %% N2E1
    ABC(G+di(30),30)     = -eta(i,j+1)/N.dx/N.dz;   % Vx at N1E1
    %% N1E1
    ABC(G+di(33),33)     = eta(i,j+1)/N.dx/N.dz;    % Vx at E1
    %% E
    ABC(G+di(32),32)     = N.Kcont/N.dz;            % P at E1
    ABC(G+di(34),34)     = eta(i,j+1)/N.dx^2;
    %% S1E1
    ABC(G+di(35),35)     = -N.Kcont/N.dz;           % P at S1E1
end
end

function [ABC] = contequat(ABC,G,di,i,j,N,B)
%% BOUNDARY COEFFIECIENTS ----------------------------------------------- %
if(i>1)&&(j>1)
    switch B.rbc
        case {0,1}
            if (j==2)&&(i==2) || (j==2)&&(i==N.nz) || (j==3)&&(i==2) || ...
                    (j==N.nx)&&(i==2) || (j==N.nx)&&(i==N.nz)
                %% ONE INTERIOR CELL ------------------------------------ %
                ABC(G-2+di(19),19)    =     N.Kbond;      % P(i,j)
                %% UPPER AND LOWER LEFT CORNER -------------------------- %
                % --------------- dP/dx=0 => P(i,j)-P(i,j+1) = 0 -------- %
                if (i==2)&&(j==2) || (i==N.nz)&&(j==2)
                    % East
                    ABC(G-2+di(34),34)=     -N.Kbond;     % P(i,j+1)
                end
                %% UPPER AND LOWER RIGHT CORNER ------------------------- %
                % --------------- dP/dx=0 => P(i,j-1)-P(i,j) = 0 -------- %
                if (i==2)&&(j==N.nx) || (i==N.nz)&&(j==N.nx)
                    % West
                    ABC(G-2+di(4),4)  =     -N.Kbond;     % P(i,j-1)
                end
            else
                %% INNER GRID POINTS ------------------------------------- %
                %% N1W1
                ABC(G-2+di(2),2)    =   -N.Kcont/N.dx;      % vx(i-1,j-1)
                ABC(G-2+di(3),3)    =   -N.Kcont/N.dz;      % vz(i-1,j-1)
                %% W1
                ABC(G-2+di(6),6)    =   N.Kcont/N.dz;       % vz(i,j-1)
                %% N1
                ABC(G-2+di(17),17)  =   N.Kcont/N.dx;       % vx(i-1,j)
            end
    end
end
end

function [di] = diagnumber(di,np,nz)


%% N1W1
di(1)       = -np*nz-np;
di(2)       = -np*nz-np+1;
di(3)       = -np*nz-np+2;
%% W1
di(4)       = -np*nz;
di(5)       = -np*nz+1;
di(6)       = -np*nz+2;
%% S1W1
di(7)       = -np*nz+np;
di(8)       = -np*nz+np+1;
di(9)       = -np*nz+np+2;
%% S2W1
di(10)      = -np*nz+2*np;
di(11)      = -np*nz+2*np+1;
di(12)      = -np*nz+2*np+2;
%% N2
di(13)      = -2*np;
di(14)      = -2*np+1;
di(15)      = -2*np+2;
%% N1
di(16)      = -np;
di(17)      = -np+1;
di(18)      = -np+2;
%% Z
di(19)      = 0;
di(20)      = 0+1;
di(21)      = 0+2;
%% S1
di(22)      = np;
di(23)      = np+1;
di(24)      = np+2;
%% S2
di(25)      = 2*np;
di(26)      = 2*np+1;
di(27)      = 2*np+2;
%% N2E
di(28)      = np*nz-2*np;
di(29)      = np*nz-2*np+1;
di(30)      = np*nz-2*np+2;
%% N1E
di(31)      = np*nz-np;
di(32)      = np*nz-np+1;
di(33)      = np*nz-np+2;
%% E
di(34)      = np*nz;
di(35)      = np*nz+1;
di(36)      = np*nz+2;
%% S1E
di(37)      = np*nz+np;
di(38)      = np*nz+np+1;
di(39)      = np*nz+np+2;
%% S2E
di(40)      = np*nz+2*np;
di(41)      = np*nz+2*np+1;
di(42)      = np*nz+2*np+2;

end

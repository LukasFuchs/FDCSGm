function [T] = Diff1D_stationary(nz,Z,k,H,rho,Ttop,Tbot)
%% Physikalischer Parameter --------------------------------------------- %
%   H       [ m ]
% Oberkruste 
%   kUC     [ W/K/m ]
%   HUC     [ W/m^3 ]
%   ucB     [ m ]
% Unterkruste 
%   kLC     [ W/K/m ]
%   HLC     [ W/m^3 ]
%   lcB     [ m ]
% Mantel   
%   kM      [ W/K/m ]
%   HM      [ W/^3 ] 
%   Ttop    [ K ]
%   Tbot    [ K ]
%% Numerische Parameter ------------------------------------------------- %
dz      =   -Z/(nz-1);
% z       =   (0:dz:-Z)';
% indUC   =   z<=-ucB; 
% indLC   =   z>-ucB & z<=-lcB; 
% indM    =   z>-lcB; 
%% Erstellung des Anfangstemperaturfeld --------------------------------- %
% rho         = zeros(nz,1);
% rho(indUC)  = rhoUC; 
% rho(indLC)  = rhoLC; 
% rho(indM)   = rhoM; 

% H           = zeros(nz,1); 
% H(indUC)    = HUC; 
% H(indLC)    = HLC; 
% H(indM)     = HM; 

% Q           = -H; 

% k           = zeros(nz,1);
% k(indUC)    = kUC; 
% k(indLC)    = kLC; 
% k(indM)     = kM; 

% Boundary conditions
T           = zeros(nz,1); 
T(1)        = Ttop; 
T(end)      = Tbot; 

%% Define Coeffizients for Matrix A and boundary conditions ------------- %
% Erstellung der durchlaufenden Indizes ----------------------------- %
Number  = zeros(nz,1);
num     = 1;
for i=1:nz
    Number(i) = num;
    num = num+1;
end

% Setup coefficient matrix A ---------------------------------------- %
a       = 1/dz^2;
% Define diagonals for matrix
diag    = zeros(nz,3);

% Inner index ------------------------------------------------------- %
iind   = Number(2:(nz-1));

diag(iind-1,1)      = a.*(k(iind) + k(iind-1))./2;
diag(iind,2)        = -a.*((k(iind+1) + k(iind))./2 + ...
                        (k(iind) + k(iind-1))./2);
diag(iind+1,3)      = a.*(k(iind+1) + k(iind))./2;

% Top and bottom index and boundary conditions ---------------------- %
bind                = Number(1,1);
diag(bind,2)        = 1;
tind                = Number(nz,1);
diag(tind,2)        = 1;

A       = spdiags(diag,[-1,0,1],nz,nz);

%% Solve system of equations -------------------------------------------- %
rhs         = T-H.*rho;

T      = A\rhs;

% q           = zeros(nz-1,1);
% zq          = zeros(nz-1,1); 
% for j=1:nz-1
%     q(j)    =   (k(j)+k(j+1))*(T(j+1)-T(j))/2/dz; 
%     zq(j)   =   dz/2*j; 
% end
% % q(1,1)      = k(1)*(T(2)-T(1))/dz;
% % q(nz,1)     = k(nz)*(T(nz)-T(nz-1))/dz;


% %% Plot solution -------------------------------------------------------- %
% % set(figure(1),'Position',[331.4,231.4,914.4,420])
% clf
% subplot(3,1,1)
% plot(T-273.15,z/1e3)
% axis ij
% 
% subplot(3,1,2)
% plot(H,z/1e3)
% axis ij
% 
% subplot(3,1,3)
% plot(q,zq/1e3)
% axis ij
% % pcolor(N.x/1e3,N.z/1e3,D.T); shading interp; colorbar
% % hold on
% % contour(N.x/1e3,N.z/1e3,D.T,100:100:1500,'k');
% % xlabel('x [km]'); ylabel('z [km]'); zlabel('Temperature [^oC]')
% % title(['Stationary temperature field for ',B.ltbc,' lateral boundary conditions'])
% % caxis([0 900])
% % set(gca,'FontWeight','Bold')

end
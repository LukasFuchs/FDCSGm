function [T] = HarroRoutine(nz,L,ucB,kUC,HUC,rhoUC,lcB,kLC,HLC,rhoLC,...
            kM,HM,rhoM,T1,T0)
% Solves the 1D steady state heat equation for 
% depth dependent heat generation rate, heat conductivity, density 
% They may be given in several layers
% for the time being they are constant in each layer, but one may
% modify the depth dependence if desired
% An exponentially decaying heat generation with depth is added to the
% heat generation rates given for each layer
% Boundary condiztion: fixed T at top and bottom

% clear all
%clf

% nz = 201;                 % no. of grid elements
nlay = 3;                 % no. of layers
ztop = [0 -ucB -lcB];    % depth of the top of each layer [m]
rho0 = [rhoUC rhoLC rhoM];  % density of layers [kg/m^3]
k0   = [kUC kLC kM];         % heat conductivity of layers [W/(K m)]
                          % Estimates from Clauser and Huenges considering
                          % composition, T and P dependence pi mal Daumen
Hbg0 = [HUC HLC HM];% const background heat generation of layers [W/kg]
% H0   = 0e-9;              % surface heat generation rate for exp depth decaying contribution
% d = 5000;                 % Decay depth of exp decaying H in [m]: H = H0*exp(-z/d)
% L = 200*1e3;              % Total depth of all layers [m]
% T1 =1282;                 % Bottom temperature
% T1 = 1315 + 273.15 + 0.5*L/1e3;
% T0 = 273.15;                   % Top temperature

ztop(nlay+1) = -L;
dz = -L/(nz-1);
z = [0:dz:-L];

% Distribute values of layers to the grid
for i = 1:nz
%     rho(i) = 0;
    k(i)   = 0;
    Hbg(i) = 0.;
    for j = 1:nlay
%         rho(i) = rho(i) + rho0(j)*(z(i)>=ztop(j))*(z(i)<ztop(j+1));
        k(i)   = k(i)   + k0(j)*(z(i)>=ztop(j))*(z(i)<ztop(j+1));
        Hbg(i) = Hbg(i) + Hbg0(j)*(z(i)>=ztop(j))*(z(i)<ztop(j+1));
    end
    if i==nz 
%         rho(i) = rho0(nlay);
        k(i) = k0(nlay);
        Hbg(i) = Hbg0(nlay);
    end
    H(i) = Hbg(i);%  + H0*exp(-z(i)/d);
end

neq = nz-2;

% Build system matrix
a = 1/2/dz^2;
b = -1/2/dz^2;

A = zeros(neq,neq);
A(1,1) = b*(k(1)+2*k(2)+k(3)); 
A(1,2) = a*(k(3)+k(2));
A(neq,neq-1) = a*(k(nz-1)+k(nz-2));
A(neq,neq) = b*(k(nz)+2*k(nz-1)+k(nz-2)); 

for j = 3:nz-2
    l = j-1;
    A(l,l) = b*(k(j+1)+2*k(j)+k(j-1));
    A(l,l+1) = a*(k(j+1)+k(j)); 
    A(l,l-1) = a*(k(j)+k(j-1));
end

% Build RHS vector
% RHS = -rho(2:nz-1).*H(2:nz-1);
RHS = -H(2:nz-1);
RHS(neq) = RHS(neq) - a*(k(nz)+k(nz-1))*T1;
RHS(1) = RHS(1) - a*(k(1)+k(2))*T0;
RHS = RHS';

% SOlve system of equatins
f = A\RHS;

% Solution back to T-field
T = zeros(nz,1);
T(1)    = T0; 
T(2:nz-1)=f;
T(nz) = T1;

% Tana    =   T0 + (rho.*H0*d^2+(qbot-rho*H0*d*exp(-h/d)+Hbg*rho*h)*z-0.5*Hbg*rho*z.^2-rho*H0*d^2*exp(-z/d))/k;

% %Plot results
% figure(1)
% subplot(3,1,1)
% plot(z/1000,T,'k','linewidth',1)
% hold on
% xlabel('Depth (km)')
% ylabel('Temperature (°C)')
% title(['T-profile for depth-dep k,H and rho, no. of layers: ' num2str(nlay)])
% 
% subplot(3,1,2)
% plot(z/1000, H,'k','linewidth',1)
% xlabel('Depth (km)')
% ylabel('Heat generation rate (W/kg)')

% q=zeros(nz,1);
% for i=2:nz-1
%     kmean=(k(i+1)+2*k(i)+k(i-1))/4;
%     q(i)=kmean*(T(i+1)-T(i-1))/2/dz;
% %     q(i)=k(i)*(T(i+1)-T(i-1))/2/dz;
% %     q(i)=k(i)*(T(i+1)-T(i))/dz;
% end
% %     q(1)= k(1)*(T(2)-T(1))/dz
%     q(1)=k(1)*(-0.5*T(3)+2*T(2)-1.5*T(1))/dz; % Looks like a factor 2 is needed. Error in Lliboutry?
% %     q(1)=k(1)*(T(4)/3-1.5*T(3)+3*T(2)-11/6*T(1))/dz;
%     q(nz)=k(nz)*(0.5*T(nz-2)-2*T(nz-1)+1.5*T(nz))/dz;
% % subplot(3,1,3)
% % plot(z/1000,q,'k','linewidth',1)
% % axis([0 L/1000 0 max(q)*1.1])
% % ylabel('Heat flux (W/m^2)')
% % xlabel('Depth (km)')

% figure(2)
% subplot(2,1,1)
% plot(z/1000,rho)
% axis([0 L/1000 0 4000])
% xlabel('Depth (km)')
% ylabel('Density (kg/m^3)')
% 
% 
% subplot(2,1,2)
% plot(z/1000,k)
% axis([0 L/1000 0 5])
% xlabel('Depth (km)')
% ylabel('Heat conductivity [W/(K m)')
% 
% disp(['heatflux top (W/m^2    =' num2str(q(1))])
% disp(['heatflux bottom (W/m^2 =' num2str(q(nz))])

end
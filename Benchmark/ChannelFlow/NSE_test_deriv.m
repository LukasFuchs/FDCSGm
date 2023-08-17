% Program to calcutlate the first derivative of the function f=k*ds/dx with
% different difference schemes. 

clear
clc
clf

% Scaling parameters
Hsc         = 400e3; 
vsc         = 1.58e-9; 
etasc       = 1e20; 

H           = 400e3/Hsc; 
z0          = 0; 
nz          = 51; 
dz          = H/(nz-1); 


z           = (z0:dz:H)';

v0          = 1.58e-9;         % Velocity at the top [m/s]
v1          = 0;                % Velocity at the bottom
v0          = v0/vsc;

eta0        = 20; 
eta1        = 25;
m           = 10^eta1/10^eta0; 

k           = logspace(eta0,eta1,nz)';
k           = k./etasc; 

eta         = 10^eta0.*exp(log(m).*z./H);
eta         = eta./etasc;

% non dimensional pressure gradient
dPdx        = -0000; 
dfdx        = ones(nz,1).*dPdx; 
dfdx(1)     = v1; 
dfdx(nz)    = v0; 

% % Exponentially increasing viscosity
% v_ana       = -v0*(m./((m-1).*exp(log(m).*z)) - 1/(m-1)) + v0;

% With horizontal pressure gradient dPdx
% v_ana       = -dPdx*H/eta0/log(m).*((z+H/(m-1)).*exp(-log(m).*z) - H/(m-1))...
%     - (H/(m-1)+H)*v0.*exp(-log(m).*z)...
%     + v0*(H/(m-1) + H);
v_ana       = -dPdx*H/(10^eta0/etasc)/log(m).*(m.^(-z./H)/(m-1).*(z.*(m-1)+H) - H/(m-1))...
    - m.^(-z./H).*m.*v0./(m-1)...
    + v0*m/(m-1);

if m == 1
%     v_ana  = v0.*z;
    v_ana = -dPdx/2/(10^eta0/etasc)*(H.*z-z.^2) + v0.*z./H;
end

% v_ana       = -1/log(m).*(1*dPdx.*m.^((-z+1)/1).*z-1*dPdx.*m.^(-z/1).*z+1^2*dPdx.*m^(-1/1)+m.^((-z+1)./1).*v0*eta0*log(m)-dPdx*1^2-v0*eta0*log(m)*m)./eta0./(m-1);

%% staggered grid
a1      = zeros(nz,1);
b1      = zeros(nz,1);
c1      = zeros(nz,1); 

for i = 2:(nz-1)
    a1(i-1) = (k(i-1)+k(i))/2/dz^2;
    b1(i)   = -(k(i-1)+2*k(i)+k(i+1))/2/dz^2;
    c1(i+1) = (k(i)+k(i+1))/2/dz^2;
end

b1(1)   = 1; 
b1(nz)  = 1; 

abc1    = [a1,b1,c1];

A1      = spdiags(abc1,[-1 0 1],nz,nz);

% spky(A1)
% rhs     = dfdx./k;

v_num   = A1\dfdx;

% STRESS CALCULATION
a12     = zeros(nz,1); 
b12     = zeros(nz,1); 
c12     = zeros(nz,1); 

for i = 2:(nz-1)
    a12(i-1)    = -(k(i-1)+k(i))/4/dz;
    b12(i)      = (k(i-1)-k(i+1))/4/dz;
    c12(i+1)    = (k(i)+k(i+1))/4/dz;
end

b12(1)      = -(k(1)+k(2))/2/dz;
b12(end)    = (k(end)+k(end-1))/2/dz;
c12(2)      = (k(1)+k(2))/2/dz;
a12(end-1)  = -(k(end)+k(end-1))/2/dz;

abc12   = [a12,b12,c12];

A12     = spdiags(abc12,[-1 0 1],nz,nz);

tauxz1  = A12*v_num;
% dtauxz1 = tauxz1(1)-tauxz1(end);
% plot(tauxz1,z)

%% Chain law
a2      = zeros(nz,1);
b2      = zeros(nz,1);
c2      = zeros(nz,1);

for i = 2:(nz-1)
    a2(i-1) = (k(i-1)+4*k(i)-k(i+1))/4/dz^2;
    b2(i)   = -2*k(i)/dz^2;
    c2(i+1) = -((k(i-1)-4*k(i)-k(i+1))/4/dz^2);
end

b2(1)   = 1;
b2(nz)  = 1;

abc2    = [a2,b2,c2];

A2      = spdiags(abc2,[-1 0 1],nz,nz);

% figure(2)
% spy(A2)

v_num2      = A2\dfdx;

%% Centered without staggered grid points
a3      = zeros(nz,1);
b3      = zeros(nz,1);
c3      = zeros(nz,1);
as      = zeros(nz,1);
bs      = zeros(nz,1); 

b3(1)   = 1; 
b3(nz)  = 1; 

as(1)   = k(1)/2/dz^2;
b3(2)   = -(k(3)+2*k(2))/4/dz^2;
c3(4)   = k(3)/4/dz^2; 

bs(nz)  = k(nz)/2/dz^2;
b3(nz-1)= -(k(nz-1)+2*k(nz))/4/dz^2;
a3(nz-3)= k(nz-2)/4/dz^2;

for i = 3:(nz-2)
    a3(i-2) = k(i-1)/4/dz^2;
    b3(i)   = -(k(i-1)+k(i+1))/4/dz^2;
    c3(i+2) =k(i+1)/4/dz^2;
end

abc3 = [a3,as,b3,bs,c3];

A3  = spdiags(abc3,[-2 -1 0 1 2],nz,nz);

% figure(2)
% spy(A3)

v_num3      = A3\dfdx;

% STRESS CALCULATION
a32         = zeros(nz,1);
b32         = zeros(nz,1); 
c32         = zeros(nz,1); 
d32         = zeros(nz,1); 
e32         = zeros(nz,1); 

for i = 2:(nz-1)
    a32(i-1) = -k(i)/2/dz;
    c32(i+1) = k(i)/2/dz;
end

b32(1)      = -3*k(1)/2/dz;
c32(2)      = 2*k(1)/dz;
d32(3)      = -k(1)/2/dz;
b32(end)    = 3*k(end)/2/dz;
a32(end-1)  = -2*k(end)/dz;
e32(end-2)  = k(end)/2/dz;

abc32       = [e32,a32,b32,c32,d32];

A32         = spdiags(abc32,[-2 -1 0 1 2],nz,nz);

% full(A32)
% spy(A32)

tauxz3      = A32*v_num3;
% plot(tauxz2,z)

% return

%% Fourth order, centered without staggered grid points
a4          = zeros(nz,1);
b4          = zeros(nz,1);
c4          = zeros(nz,1);
d4          = zeros(nz,1);
e4          = zeros(nz,1);
f4          = zeros(nz,1);
g4          = zeros(nz,1);
h4          = zeros(nz,1);
j4          = zeros(nz,1);
m4          = zeros(nz,1);
z4          = zeros(nz,1);

e4(1)       = 1;
e4(nz)      = 1;

a           = 1/12;
b           = -2/3;

ap          = -25/12;
bp          = 4;
cp          = -3;
dp          = 4/3;
ep          = -1/4;

a2p         = -1/4;
b2p         = -5/6;
c2p         = 3/2;
d2p         = -1/2;
e2p         = 1/12;


for i = 5:(nz-4)
    a4(i-4) = a^2*k(i-2);
    b4(i-3) = a*b*(k(i-2) + k(i-1));
    c4(i-2) = b^2*k(i-1);
    d4(i-1) = -a*b*(k(i-2) + k(i+1));
    e4(i)   = -(a^2*(k(i-2) + k(i+2)) + b^2*(k(i-1) + k(i+1)));
    f4(i+1) = -a*b*(k(i-1) + k(i+2));
    g4(i+2) = b^2*k(i+1);
    h4(i+3) = a*b*(k(i+1) + k(i+2));
    j4(i+4) = a^2*k(i+2);
end

% Boundaries
% for i = 4
b4(1)   = a*(a2p*k(2) + b*k(3));
c4(2)   = a*b2p*k(2) + b^2*k(3);
d4(3)   = a*(k(2)*c2p - b*k(5));
e4(4)   = a*(k(2)*d2p - a*k(6)) - b^2*(k(3)+k(5));
f4(5)   = a*(k(2)*e2p - b*(k(3) + k(6)));
g4(6)   = b^2*k(5);
h4(7)   = a*b*(k(5) + k(6));
j4(8)   = a^2*k(6);
% end

% for i = nz-3
a4(nz-7)    = a^2*k(nz-5);
b4(nz-6)    = a*b*(k(nz-5) + k(nz-4));
c4(nz-5)    = b^2*k(nz-4);
d4(nz-4)    = -a*(b*(k(nz-5) + k(nz-2)) - e2p*k(nz-1));
e4(nz-3)    = -(a*(a*k(nz-5) - d2p*k(nz-1)) + b^2*(k(nz-4) + k(nz-2)));
f4(nz-2)    = a*(c2p*k(nz-1) - b*k(nz-4));
g4(nz-1)    = b^2*k(nz-2) + a*b2p*k(nz-1);
h4(nz)      = a*(b*k(nz-2) + a2p*k(nz-1));
% end

% for i = 3
c4(1)   = a*ap*k(1) + b*a2p*k(2);
d4(2)   = a*(bp*k(1) - b*k(4)) + b*b2p*k(2);
e4(3)   = a*(cp*k(1) - a*k(5)) + b*(c2p*k(2) - b*k(4));
f4(4)   = a*(dp*k(1) - b*k(5)) + b*d2p*k(2);
g4(5)   = a*ep*k(1) + b*(e2p*k(2) + b*k(4));
h4(6)   = a*b*(k(4) + k(5));
j4(7)   = a^2*k(5);
% end

% for i = nz-2
a4(nz-6) = a^2*k(nz-4);
b4(nz-5) = a*b*(k(nz-4) + k(nz-3));
c4(nz-4) = b*(b*k(nz-3) + e2p*k(nz-1)) + a*ep*k(nz);
d4(nz-3) = -(b*(a*k(nz-4) - d2p*k(nz-1)) - a*dp*k(nz));
e4(nz-2) = -(a*(a*k(nz-4) - cp*k(nz)) + b*(b*k(nz-3) - c2p*k(nz-1)));
f4(nz-1) = -(a*(b*k(nz-3) - bp*k(nz)) - b*b2p*k(nz-1));
g4(nz)   = b*a2p*k(nz-1) + a*ap*k(nz);
% end

% for i = 2
d4(1) = a2p*(ap*k(1) + b2p*k(2)) + a*c2p*k(3);
e4(2) = bp*a2p*k(1) + b2p^2*k(2) + b*c2p*k(3) + a*d2p*k(4) ;
f4(3) = cp*a2p*k(1) + b2p*c2p*k(2) + b*d2p*k(4) + a*e2p*k(5);
g4(4) = dp*a2p*k(1) + b2p*d2p*k(2) - b*(c2p*k(3) - e2p*k(5));
h4(5) = ep*a2p*k(1) + b2p*e2p*k(2) - a*c2p*k(3) - b*d2p*k(4);
j4(6) = -(a*d2p*k(4) + b*e2p*k(5));
m4(7) = -a*e2p*k(5);
% end

% for i = nz-1
z4(nz-6) = -(a*e2p*k(nz-4));
a4(nz-5) = -(b*e2p*k(nz-4) + a*d2p*k(nz-3));
b4(nz-4) = -(b*d2p*k(nz-3) + a*c2p*k(nz-2) - b2p*e2p*k(nz-1) - ep*a2p*k(nz));
c4(nz-3) = b*(e2p*k(nz-4) - c2p*k(nz-2)) + b2p*d2p*k(nz-1) + dp*a2p*k(nz);
d4(nz-2) = a*e2p*k(nz-4) + b*d2p*k(nz-3) + c2p*b2p*k(nz-1) + cp*a2p*k(nz);
e4(nz-1) = a*d2p*k(nz-3) + b*c2p*k(nz-2) + b2p^2*k(nz-1) + bp*a2p*k(nz);
f4(nz)   = a*c2p*k(nz-2) + a2p*b2p*k(nz-1) + ap*a2p*k(nz);
% end

z4      = z4./dz^2;
a4      = a4./dz^2;
b4      = b4./dz^2;
c4      = c4./dz^2;
d4      = d4./dz^2;
e4      = e4./dz^2;
f4      = f4./dz^2;
g4      = g4./dz^2;
h4      = h4./dz^2;
j4      = j4./dz^2;
m4      = m4./dz^2;

e4(1)   = 1; 
e4(nz)  = 1;

abc4    = [z4,a4,b4,c4,d4,e4,f4,g4,h4,j4,m4];

A4      = spdiags(abc4,[-5 -4 -3 -2 -1 0 1 2 3 4 5],nz,nz);

% figure(2)
% spy(A4)

v_num4      = A4\dfdx;

% STRESS CALCULATION
x42         = zeros(nz,1);
z42         = zeros(nz,1);
a42         = zeros(nz,1);
b42         = zeros(nz,1);
c42         = zeros(nz,1);
d42         = zeros(nz,1);
e42         = zeros(nz,1);
f42         = zeros(nz,1);
g42         = zeros(nz,1);

for i = 3:(nz-2)
    a42(i-2)    = k(i)/12/dz;
    b42(i-1)    = -2*k(i)/3/dz;
    d42(i+1)    = -b42(i-1);
    e42(i+2)    = -a42(i-2);
end
for i = 2
    b42(i-1)    = -k(i)/4/dz;
    c42(i)      = -5*k(i)/6/dz;
    d42(i+1)    = 3*k(i)/2/dz;
    e42(i+2)    = -k(i)/2/dz;
    f42(i+3)    = k(i)/12/dz;
end
for i = (nz-1)
   d42(i+1)     = k(i)/4/dz;
   c42(i)       = 5*k(i)/6/dz;
   b42(i-1)     = -3*k(i)/2/dz;
   a42(i-2)     = k(i)/2/dz;
   z42(i-3)     = -k(i)/12/dz; 
end
for i = 1
    c42(i)      = -25*k(i)/12/dz;
    d42(i+1)    = 4*k(i)/dz;
    e42(i+2)    = -3*k(i)/dz;
    f42(i+3)    = 4*k(i)/3/dz;
    g42(i+4)    = -k(i)/4/dz;
end
for i = nz
    c42(i)      = 25*k(i)/12/dz;
    b42(i-1)    = -4*k(i)/dz;
    a42(i-2)    = 3*k(i)/dz;
    z42(i-3)    = -4*k(i)/3/dz;
    x42(i-4)    = k(i)/4/dz;
end

abc42       = [x42,z42,a42,b42,c42,d42,e42,f42,g42];
A42         = spdiags(abc42,[-4 -3 -2 -1 0 1 2 3 4],nz,nz); 
% full(A42)
% spy(A42)

tauxz4      = A42*v_num4;
% plot(tauxz1,z,'k',tauxz4,z,'m')
% return


% v_num       = vnum./v0;
% v_num2      = v_num2./v0;
% v_num3      = v_num3./v0;
% v_num4      = v_num4./v0;
% v_ana       = v_ana./v0;

d_num1      = sqrt((v_num-v_ana).^2./max(v_ana).^2);
d_num2      = sqrt((v_num2-v_ana).^2./max(v_ana).^2);
d_num3      = sqrt((v_num3-v_ana).^2./max(v_ana).^2);
d_num4      = sqrt((v_num4-v_ana).^2./max(v_ana).^2);


%%
nmsize      = 5; 
% z           = z./1000;
figure(1)
subplot(3,2,1)
% plot(z,s,'b',z,dsdx,'r',z,k,'g',z,v_num,'k',z,v_num2,'--k',z,v_num3,'k*')
plot(v_ana,z,'k-',v_num,z,'k*',v_num2,z,'rd',v_num3,z,'gs',v_num4,z,'mo','MarkerSize',nmsize)
grid on
% axis([0 max(v_ana) 0 H])
% legend('s','dsdx','k','v_{num}','v_{num2}','v_{num3}','Location','Best')
% legend('v_{ana}','v_{steg}','v_{Chain}','v_{centr2nd}','v_{centr4th}','Location','Best')
ylabel('z')
xlabel('\itv')
title('Velocity [non dim]')

subplot(3,2,3) 
semilogx(k,z,'b',eta,z,'d','MarkerSize',nmsize)
ylabel('z')
xlabel('\it\eta')
title('Viscosity')
grid on

subplot(1,2,2)
semilogy(z,d_num1,'k*',z,d_num2,'rd',z,d_num3,'gs',z,d_num4,'mo','MarkerSize',nmsize)
grid on
legend('\Delta_{stag}','\Delta_{Chain}','\Delta_{centr2nd}','\Delta_{centr4th}','Location','Best')
xlabel('z')
ylabel('\Delta')
title('Square root error')

% figure(2)
subplot(3,2,5)
plot(tauxz1,z,'k',tauxz4,z,'m--')
grid on
xlabel('\tau_{xz}')
ylabel('z')
title('Shear stress')
legend('stag','centr4th','Location','Best')





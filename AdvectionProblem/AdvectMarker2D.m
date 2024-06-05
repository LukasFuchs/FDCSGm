function [XM,ZM]   = AdvectMarker2D(M,XM,ZM,dt,vx,vz)

% x = dt*v(n-2) --------------------------------------------------------- %
vxm = interp2(M.x,M.z,vx,XM,ZM,'linear');
vzm = interp2(M.x,M.z,vz,XM,ZM,'linear');
x1  = dt.*vxm;
z1  = dt.*vzm;

% x = dt*v(n-1)  -------------------------------------------------------- %
X   = XM + x1./2;
Z   = ZM + z1./2;
X(X>M.L)    = X(X>M.L) - (max(X(X>M.L))-M.L);
X(X<0)      = X(X<0) - (min(X(X<0))-0);
Z(Z>0)      = Z(Z>0) - (max(Z(Z>0))-0);
Z(Z<M.H)    = Z(Z<M.H) - (min(Z(Z<M.H))-M.H);
vxm = interp2(M.x,M.z,vx,X,Z,'linear');
vzm = interp2(M.x,M.z,vz,X,Z,'linear');
x2  = dt.*vxm;
z2  = dt.*vzm;

% x = dt*v(n-1)  -------------------------------------------------------- %
X   = XM + x2./2;
Z   = ZM + z2./2;
X(X>M.L)    = X(X>M.L) - (max(X(X>M.L))-M.L);
X(X<0)      = X(X<0) - (min(X(X<0))-0);
Z(Z>0)      = Z(Z>0) - (max(Z(Z>0))-0);
Z(Z<M.H)    = Z(Z<M.H) - (min(Z(Z<M.H))-M.H);
vxm = interp2(M.x,M.z,vx,X,Z,'linear');
vzm = interp2(M.x,M.z,vz,X,Z,'linear');
x3  = dt.*vxm;
z3  = dt.*vzm;

% x = dt*v(n)  ---------------------------------------------------------- %
X   = XM + x3;
Z   = ZM + z3;
X(X>M.L)    = X(X>M.L) - (max(X(X>M.L))-M.L);
X(X<0)      = X(X<0) - (min(X(X<0))-0);
Z(Z>0)      = Z(Z>0) - (max(Z(Z>0))-0);
Z(Z<M.H)    = Z(Z<M.H) - (min(Z(Z<M.H))-M.H);
vxm = interp2(M.x,M.z,vx,X,Z,'linear');
vzm = interp2(M.x,M.z,vz,X,Z,'linear');
x4  = dt.*vxm;
z4  = dt.*vzm;

XM   = XM + (x1+2.*(x2+x3)+x4)./6;
ZM   = ZM + (z1+2.*(z2+z3)+z4)./6;

XM(XM>M.L)  = 0 + (XM(XM>M.L)-M.L);
XM(XM<0)    = M.L + (XM(XM<0)-0);
ZM(ZM>0)    = M.H + (ZM(ZM>0)-0);
ZM(ZM<M.H)  = 0 + (ZM(ZM<M.H)-M.H);
end
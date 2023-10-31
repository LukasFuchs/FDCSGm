function ID = VelGradTensor(ID,N)

% dv_x/dx --------------------------------------------------------------- %
ID.dvxdx    =   velgrad(ID.vx,N.dx,'hor');
% dv_z/dz --------------------------------------------------------------- %
ID.dvzdz    =   velgrad(ID.vz,N.dz,'vert');
% dv_x/dz --------------------------------------------------------------- %
ID.dvxdz    =   velgrad(ID.vx,N.dz,'vert');
% dv_z/dx --------------------------------------------------------------- %
ID.dvzdx    =   velgrad(ID.vz,N.dx,'hor');

% Divergence ------------------------------------------------------------ %
ID.div      =   ID.dvxdx+ID.dvzdz;
% Vorticity ------------------------------------------------------------- %
ID.rot      =   ID.dvxdz-ID.dvzdx;


end
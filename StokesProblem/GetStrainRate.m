function [ID] = GetStrainRate(ID,N)

% Invarianten der Dehnungsrate und der Spannung ------------------------- %
ID.dvxdx    =   velgrad(ID.vx,N.dx,'hor');
ID.dvzdz    =   velgrad(ID.vz,N.dz,'vert');
ID.dvzdx    =   velgrad(ID.vz,N.dx,'hor');
ID.dvxdz    =   velgrad(ID.vx,N.dz,'vert');

ID.eII     	=   sqrt(0.5*(ID.dvxdx.^2 + ID.dvzdz.^2) + ...
                    0.25.*(ID.dvxdz + ID.dvzdx).^2);


% Strain rate tensor ---------------------------------------------------- %
% e_ij = 0.5 * ( dv_i/dx_j + dv_j/dx_i )
% ----------------------------------------------------------------------- %

% keyboard

end
function [D] = GetStrainRateStress(D,N,M)

% Invarianten der Dehnungsrate und der Spannung ------------------------- %
% ID.dvxdx    =   velgrad(ID.vx,N.dx,'hor');
% ID.dvzdz    =   velgrad(ID.vz,N.dz,'vert');
% ID.dvzdx    =   velgrad(ID.vz,N.dx,'hor');
% ID.dvxdz    =   velgrad(ID.vx,N.dz,'vert');
% 
% ID.eII     	=   sqrt(0.5*(ID.dvxdx.^2 + ID.dvzdz.^2) + ...
%                     0.25.*(ID.dvxdz + ID.dvzdx).^2);


% Strain rate tensor ---------------------------------------------------- %
% e_ij = 0.5 * ( dv_i/dx_j + dv_j/dx_i )
% ----------------------------------------------------------------------- %
iind    =   1:N.nz1;
jind    =   1:N.nx1; 

D.exx(iind,jind)    = (D.vx(iind+1,jind)-D.vx(iind,jind))./N.dx;
D.ezz(iind,jind)    = (D.vz(iind+1,jind)-D.vz(iind,jind))./N.dz;
D.exz(iind,jind)    = 

keyboard

end
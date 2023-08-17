function [ID] = InterpStaggered(D,ID,N,param)
% Function to interpolate a property from the staggered grid on the regular
% grid
% ======================================================================= %

switch lower(param)
    case 'velocity'
        % Interal grid points ------------------------------------------- %
        ID.vx(ID.indi,ID.indj)    = ...
            (D.vx(ID.indi-1,ID.indj) + D.vx(ID.indi,ID.indj)) / 2;
        ID.vz(ID.indi,ID.indj)    = ...
            (D.vz(ID.indi,ID.indj-1) + D.vz(ID.indi,ID.indj)) / 2;
        
        % Boundary nodes ------------------------------------------------ %
        ID.vx(1,:)          = D.vx(1,:) + (D.vx(1,:) - D.vx(2,:))/2;
        ID.vx(ID.indi,1)    = (D.vx(ID.indi-1,1)+D.vx(ID.indi,1))/2;
        ID.vx(ID.indi,N.nx) = (D.vx(ID.indi-1,N.nx)+D.vx(ID.indi,N.nx))/2;
        ID.vx(N.nz,:)       = D.vx(N.nz1,:) + (D.vx(N.nz1,:) - D.vx(N.nz1-1,:))/2;        
        
        ID.vz(:,1)          = D.vz(:,1) + (D.vz(:,1) - D.vz(:,2))/2;
        ID.vz(1,ID.indj)    = (D.vz(1,ID.indj-1)+D.vz(1,ID.indj))/2;
        ID.vz(N.nz,ID.indj) = (D.vz(N.nz,ID.indj-1)+D.vz(N.nz,ID.indj))/2;
        ID.vz(:,N.nx)       = D.vz(:,N.nx1) + (D.vz(:,N.nx1) - D.vz(:,(N.nx1-1))) / 2;       

        ID.v            =   sqrt(ID.vx.^2 + ID.vz.^2);
end

end
function [dvd] = velgrad(v,d,param)
% Function to calculate the velocity gradients on the regular nodal points.
% Velocity gradients are used to calculate the second invariant of the
% strain rate, which is used for the composite viscosity.
%
% Ver. 1.0. - 13.08.14
% ----------------------------------------------------------------------- %

nz      =   size(v,1);
nx      =   size(v,2);

dvd     =   zeros(nz,nx);

switch lower(param)
    
    % HORIZONTAL DERIVATIVE
    case 'hor'
        indi    =   1:nz; 
        indj    =   2:(nx-1); 

        % Internal grid points ------------------------------------------ %
        dvd(indi,indj)  = (v(indi,indj+1)-v(indi,indj-1))/2/d;
        
        % Boundary nodes ------------------------------------------------ %
        dvd(:,1)    = (-1.5.*v(:,1) + 2.*v(:,2) - 0.5.*v(:,3))./d;
        dvd(:,nx)   = (1.5.*v(:,nx) - 2.*v(:,(nx-1)) + 0.5.*v(:,(nx-2)))./d;
        
        % VERTICAL DERIVATIVE
    case 'vert'
        indi    =   2:(nz-1); 
        indj    =   1:nx; 
        
        % Internal grid points ------------------------------------------ %
        dvd(indi,indj)    = (v(indi+1,indj)-v(indi-1,indj))/2/d;
        
        % Boundary nodes ------------------------------------------------ %
        dvd(1,:)    = (-1.5.*v(1,:) + 2.*v(2,:) - 0.5.*v(3,:))/d;
        dvd(nz,:)   = (1.5.*v(nz,:) - 2.*v((nz-1),:) + 0.5.*v((nz-2),:))/d;
end


end
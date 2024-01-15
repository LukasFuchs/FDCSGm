function [D] = RHEOLOGY(D,G,R,T,GC,N,M)
% ======================================================================= %
% T     - Temperature                   [ K ]                             %
% eII   - Total strain rate             [ s^-1 ]                          %
% r     - grain size                    [ m ]                             %
%                                                                         %
% ======================================================================= %

D.tauII0    =   0;
D.epsl      =   D.eII;

for j = 1:R.nrheo
    % Rheological iteration --------------------------------------------- %
    % Get effective viscosity ------------------------------------------- %
    % Dislocation creep viscosity --------------------------------------- %
    D           =   DislocationCreep(D,R,T,GC);
    
    % Diffusion creep viscosity ----------------------------------------- %
    D           =   DiffusionCreep(D,R,T,GC);
    
    % Effective Viscosity ----------------------------------------------- %
    D           =   EffectiveViscosity(D,R);
    
    % Get stress -------------------------------------------------------- %
    D.tauII     =   2.*D.eta.*D.eII;            % [ Pa ]
    D.dtau      =   D.tauII0 - D.tauII;
    D.tauII0    =   D.tauII;
    
    D.epsln     =   D.tauII./2./D.etal;         % Disl. strain rate
    D.epsf      =   D.tauII./2./D.etaf;         % Diff. strain rate
    D.epstot    =   D.epsln + D.epsf;           % Total strain rate
    D.deps      =   D.epsl - D.epsln;
    
    % Steady State Grain Size ------------------------------------------- %
    switch lower(R.Grains)
        case 'steadystate'
            D   =   SteadyStateGrainSize(D,R,G,GC,T);
    end
    
    switch lower(R.plast)
        case {'drucker-prager','yes'}
            D.epstot    =   max(D.epstot,D.eII);
    end
    
    D.epsl  =   D.epsl*R.itfac + (1.0-R.itfac)*D.epsln;
    
    if(N.debug==1)
        figure(100)
        clf
        subplot(2,2,1)
        semilogx(D.d,M.z./1e3,'b-')
        hold on; box on; axis square
        xlabel('R [ m ]'); ylabel('Depth [ km ]')
        subplot(2,2,2)
        semilogx(D.eta,M.z./1e3,'k-',D.etal,M.z./1e3,...
            'r--',D.etaf,M.z./1e3,'b--')
        xlabel('\eta [ Pa s ]'); ylabel('Depth [km]')
        %         axis([M.H 0 1e15 1e30])
        box on; axis square
        subplot(2,2,3)
        % semilogx(abs(D.dtau),M.z./1e3,'b-'); box on
        plot(D.tauII./1e6,M.z./1e3)
        axis square
        xlabel('\tau [ MPa ]'); ylabel('Depth [ km ]')
        subplot(2,2,4)
        %  semilogx(D.deps,M.z./1e3,'r-'); box on
        semilogx(D.epstot,M.z./1e3,'k-',D.epsl,M.z./1e3,...
            'r--',D.epsf,M.z./1e3,'b--')
        axis square
        xlabel('\epsilon'); ylabel('Depth [ km ]')
%         pause(.5)
    end
    
end

end

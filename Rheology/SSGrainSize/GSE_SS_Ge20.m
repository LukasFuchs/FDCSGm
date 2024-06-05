function [D] = GSE_SS_Ge20(D,G,GC,T)
% All constants are used in SI units
%                                                                         %
% ======================================================================= %

D.lab       = 'Ge20';

D.d         =   D.d.*sqrt(pi/2);        % Initial interface roughness

% Interface roughness growth -------------------------------------------- %
GI      =   G.Ge20{5}.*exp(-G.Ge20{2}./(GC.RG.*T.T));   % [ m^p s^-1 ]
% ----------------------------------------------------------------------- %

% Constants ------------------------------------------------------------- %
f0      =   0.001;
fI      =   f0*exp(-2.*(T.T./1000).^2.9);
DI      =   fI./G.Ge20{7}./G.Ge20{6};
CI      =   G.Ge20{6}.*GI./G.Ge20{4};               % [ m^p s^-1 ]
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% % Assuming both mechanisms act simultaneously---------------------------- %
% for i=1:100
%     % Dislocation creep viscosity --------------------------------------- %
%     D.etal = 0.5*R.Al^(-1/R.n).*...
%         exp(R.Ql./GC.RG./T./R.n).*D.epsss.^(1/R.n-1);       % [ Pa s ]
%     
%     % Diffusion creep viscosity ----------------------------------------- %
%     % Grain size needs to be in m! -------------------------------------- %
%     D.etaf = 0.5*R.Af^(-1).*...
%         exp(R.Qf./GC.RG./T).*D.Rinf.^R.m;                   % [ Pa s ]
%     
%     % Effective viscosity ----------------------------------------------- %
%     D.etass     = (1./D.etaf + 1./D.etal).^(-1);
%     
%     % Get stress -------------------------------------------------------- %
%     D.tauss     = 2.*D.etass.*E;                            % [ Pa ]
    
    psi         = 2.*D.tauII.*D.epsln;
    
    % Grain size -------------------------------------------------------- %           
    rinf        = (CI./DI).^(1/(G.Ge20{4}+1)).*...
        (1./psi).^(1/(G.Ge20{4}+1));
    
    D.d         =   D.d.*G.itfac + rinf.*(1-G.itfac);
%     D.rinf      = D.rinf.*0.9 + rinf.*(1-0.9);
    
    D.d         =   D.d./sqrt(pi/2);
%     D.Rinf      = D.rinf./sqrt(pi/2);
    
%     D.efss      = D.tauss./2./D.etaf;
%     D.epsss     = D.tauss./2./D.etal;
%     D.etotss    = D.efss + D.epsss;
    
%     switch N.debug
%         case 3
%             if (max(size(D.Rinf))>1)
%                 error('ERROR! Dimension for debug check to large to plot!')
%             end
%             figure(1)
%             subplot(2,2,1)
%             plot(i,log10(D.etass(:)),'k+',...
%                 i,log10(D.etaf(:)),'b.',...
%                 i,log10(D.etal(:)),'ro')
%             title('Viscosity')
%             hold on
%             subplot(2,2,2)
%             plot(i,log10(D.tauss(:)),'k.')
%             title('Stress')
%             hold on
%             subplot(2,2,3)
%             plot(i,log10(D.efss(:)),'b.',...
%                 i,log10(D.epsss(:)),'ro',...
%                 i,log10(D.etotss(:)),'k+')
%             title('Strain Rate')
%             hold on
%             subplot(2,2,4)
%             semilogy(i,D.rinf(:),'k.',i,D.Rinf(:),'ro')
%             title('Grain Size')
%             hold on           
%             pause(.1)
%     end
% end
% ======================================================================= %
% ============================ END ====================================== %
% ======================================================================= %
end
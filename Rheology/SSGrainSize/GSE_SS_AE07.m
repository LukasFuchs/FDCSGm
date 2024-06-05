function [D] = GSE_SS_AE07(D,G,GC,R,T)
% ======================================================================= %
% All constants are used in SI units                                      %
%                                                                         %
% ======================================================================= %

switch lower(G.GSE)
    case 'ae07'
        D.lab       ='AE07';
        % Grain growth -------------------------------------------------- %
        k       = G.AE07{1}.*exp(-G.AE07{2}./(GC.RG.*T.T));   % [ m^p s^-1 ]
        
        % Kinetic term for grain coarsening ----------------------------- %
        p       = G.AE07{3};                % Grain growth exponent
        
        % Constants for grain size reduction ---------------------------- %
        gamma   = G.AE07{4};                                % [ J m^-2 ]
        lambda  = G.AE07{5};
        c       = G.AE07{6};
    case 'be09'
        D.lab       ='Be09';
        % Grain growth -------------------------------------------------- %
        k       = G.Be09{1}.*exp(-G.Be09{2}./(GC.RG.*T.T));   % [ m^p s^-1 ]
        
        % Kinetic term for grain coarsening ----------------------------- %
        p       = G.Be09{3};                % Grain growth exponent
        
        % Constants for grain size reduction ---------------------------- %
        gamma   = G.Be09{4};                                % [ J m^-2 ]
        lambda  = G.Be09{5};
        c       = G.Be09{6};
end
% ----------------------------------------------------------------------- %

% Constants ------------------------------------------------------------- %
% The use of SI units in the rheology, requires tau to be in Pa.
% This leeds to a grain size in m and constants are scaled accordingly.
C1      = 2*lambda/c/gamma;                 % [ m^2/J ]
A1      = R.Al.*exp(-R.Ql./GC.RG./T.T);       % [ s^-1 Pa^-n ]
X1      = k.*A1.^(1./R.n)./p./C1;            % [ m^(p+1) / s^((n+1)/n) ]
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% Assuming both mechanisms act simultaneously---------------------------- %
% for i=1:200
%     % Dislocation creep viscosity --------------------------------------- %
%     D.etal    = 0.5*R.Al^(-1/R.n).*...
%         exp(R.Ql./GC.RG./T./R.n).*D.epsss.^(1/R.n-1);       % [ Pa s ]
%     
%     % Diffusion creep viscosity ----------------------------------------- %
%     % Grain size needs to be in m! -------------------------------------- %
%     D.etaf    = 0.5*R.Af^(-1).*...
%         exp(R.Qf./GC.RG./T).*D.Rinf.^R.m;                   % [ Pa s ]
%     
%     % Effective viscosity ----------------------------------------------- %
%     D.etass     = (1./D.etaf + 1./D.etal).^(-1);
    
%     % Get stress -------------------------------------------------------- %
%     D.tauss     = 2.*D.etass.*E;                            % [ Pa ]
    
    % Grain size -------------------------------------------------------- %
    Rinf        = X1.^(1./(p+1)).*...
        D.epsln.^(-(R.n+1)./(R.n.*(p+1)));                    % [ m ]
    
    D.d      =  D.d.*G.itfac + Rinf.*(1-G.itfac);
    
%     D.efss      = D.tauss./2./D.etaf;
%     D.epsss     = D.tauss./2./D.etal;
%     D.etotss    = D.efss + D.epsss;
    
%     switch N.debug
%         case 3
%             if (max(size(D.Rinf))>1)
%                 error('ERROR! Dimension for debug check to large to plot!')
%             end
%             figure(101)
%             subplot(2,2,1)
%             plot(i,log10(D.etass(:)),'k+',...
%                 i,log10(etaf(:)),'b.',...
%                 i,log10(etal(:)),'ro')
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
%             semilogy(i,D.Rinf(:),'k.',i,Rinf,'ro')
%             title('Grain Size')
%             hold on
%             pause(.1)
%     end
% end
% ======================================================================= %
% ============================ END ====================================== %
% ======================================================================= %
end
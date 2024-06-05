function [D] = GSE_SS_DAN17(D,G,GC,R,T)
% ======================================================================= %
% All constants are used in SI units                                      %
%                                                                         %
% ======================================================================= %

D.lab       = 'Dan17';

% Grain growth ---------------------------------------------------------- %
KG      = G.Dan17{1}.*exp(-G.Dan17{2}./(GC.RG.*T.T));         % [ m^p s^-1 ]
% ----------------------------------------------------------------------- %

% Constants ------------------------------------------------------------- %
C1      = 2*G.Dan17{5}/G.Dan17{6}/G.Dan17{4};   % [ m^2/J ]
A1      = R.Al.*exp(-R.Ql./GC.RG./T.T);         % [ s^-1 Pa^-n ]
X1      = KG.*A1.^(1./R.n)./G.Dan17{3}./C1;     % [ m^(p+1) / s^((n+1)/n) ]
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% Assuming both mechanisms act simultaneously---------------------------- %
% for i=1:200
%     % Dislocation creep viscosity --------------------------------------- %
%     etal    = 0.5*R.Al^(-1/R.n).*...
%         exp(R.Ql./GC.RG./T./R.n).*D.epsss.^(1/R.n-1);       % [ Pa s ]
%     
%     % Diffusion creep viscosity ----------------------------------------- %
%     % Grain size needs to be in m! -------------------------------------- %
%     etaf    = 0.5*R.Af^(-1).*...
%         exp(R.Qf./GC.RG./T).*D.Rinf.^R.m;                   % [ Pa s ]
%     
%     % Effective viscosity ----------------------------------------------- %
%     D.etass     = (1./etaf + 1./etal).^(-1);
%     
%     % Get stress -------------------------------------------------------- %
%     D.tauss     = 2.*D.etass.*E;                            % [ Pa ]
    
    % Grain size -------------------------------------------------------- %
    Rinf        = X1.^(1/(G.Dan17{3}+1)).*...
        D.epsln.^(-(R.n+1)./(R.n.*(G.Dan17{3}+1)));           % [ m ]
    
    % Rinf        = X2.^(1/(G.Dan17{3}+1)).*D.tauss^((R.n+1)/(G.Dan17{3}+1));   % [ m ]
    
    D.d      =  D.d.*G.itfac + Rinf.*(1-G.itfac);
    
%     D.efss      = D.tauss./2./etaf;
%     D.epsss     = D.tauss./2./etal;
%     D.etotss    = D.efss + D.epsss;
    
%     switch N.debug
%         case 3
%             if (max(size(D.Rinf))>1)
%                 error('ERROR! Dimension for debug check to large to plot!')
%             end
%             figure(1)
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
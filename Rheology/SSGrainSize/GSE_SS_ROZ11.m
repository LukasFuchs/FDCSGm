function [D] = GSE_SS_ROZ11(D,G,GC,R,T)
% ======================================================================= %
% All constants are used in SI units                                      %
%                                                                         %
% ======================================================================= %

D.lab       ='Ro11';

% Grain growth ---------------------------------------------------------- %
KG      =   G.Ro11{1}.*exp(-G.Ro11{2}./(GC.RG.*T.T));   % [ m^p s^-1 ]
% ----------------------------------------------------------------------- %

% partitioning parameter f ---------------------------------------------- %
f       =   exp(-G.Ro11{7}.*(T.T./1000).^G.Ro11{8});
% ----------------------------------------------------------------------- %

% Constants ------------------------------------------------------------- %
C1      =   2.*f.*G.Ro11{6}./3./G.Ro11{4}./G.Ro11{5};   % [ m^2/J ]
A1      =   R.Al.*exp(-R.Ql./GC.RG./T.T);               % [ s^-1 Pa^-n ]
X       =   (KG.*A1.^(1./R.n)./G.Ro11{3}./C1);          % [ m^(p+1) / s^((n+1)/n) ]
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% % Assuming both mechanisms act simultaneously---------------------------- %
% for i=1:200
%     % Dislocation creep viscosity --------------------------------------- %
%     etal = 0.5*R.Al^(-1/R.n).*...
%         exp(R.Ql./GC.RG./T./R.n).*D.epsss.^(1/R.n-1);       % [ Pa s ]
%     
%     % Diffusion creep viscosity ----------------------------------------- %
%     % Grain size needs to be in m! -------------------------------------- %
%     etaf = 0.5*R.Af^(-1).*...
%         exp(R.Qf./GC.RG./T).*D.Rinf.^R.m;                   % [ Pa s ]
%     
%     % Effective viscosity ----------------------------------------------- %
%     D.etass     = (1./etaf + 1./etal).^(-1);
%     
%     % Get stress -------------------------------------------------------- %
%     D.tauss     = 2.*D.etass.*E;                            % [ Pa ]
    
    % Grain size -------------------------------------------------------- %
    Rinf        = X.^(1./(G.Ro11{3}+1)).*...
        D.epsln.^(-(R.n+1)./(R.n.*(G.Ro11{3}+1)));          % [ m ]
    
    %     Rinf        = ((3*G.Ro11{4}.*G2.*G.Ro11{5})/(G.Ro11{3}.*f.*A2.*G.Ro11{6})).^(1/(G.Ro11{3}+1)).*...
    %         (D.tauss./1e6).^(-((R.n+1)/(G.Ro11{3}+1)))./1e6;          % [ m ]
    
    D.d      =  D.d.*G.itfac + Rinf.*(1-G.itfac);
%     D.Rinf      = D.Rinf.*0.9 + Rinf.*(1-0.9);
    
%     D.efss      = D.tauss./2./etaf;
%     D.epsss     = D.tauss./2./etal;
%     D.etotss    = D.efss + D.epsss;
%     
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
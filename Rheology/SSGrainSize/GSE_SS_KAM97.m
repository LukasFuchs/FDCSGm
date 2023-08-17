function [D] = GSE_SS_KAM97(D,GC,G,T)
% ======================================================================= %
%                                                                         %
% ======================================================================= %

D.lab       = 'Kam97';

% Grain growth ---------------------------------------------------------- %
k           = G.Kam97{1}.*exp(-G.Kam97{2}/GC.RG./T.T);      % [ m^2/s ]
% ----------------------------------------------------------------------- %

% Scale for strain rate and stress for second invariant use ------------- %
esc         = 2/sqrt(3);
ssc         = sqrt(3);
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% Assuming both mechanisms act simultaneously
% for i=1:200
%     % Dislocation creep viscosity --------------------------------------- %
%     etal = 0.5*R.Al^(-1/R.n)*...
%         exp(R.Ql/GC.RG./T./R.n).*D.epsss.^(1/R.n-1);    % [ Pa s ]
%     
%     % Diffusion creep viscosity ----------------------------------------- %
%     % Grain size needs to be in m! -------------------------------------- %
%     etaf = 0.5*R.Af^(-1)*...
%         exp(R.Qf/GC.RG./T).*(D.Rinf).^R.m;              % [ Pa s ]
%     
%     % Effective viscosity ----------------------------------------------- %
%     D.etass     = (1./etaf + 1./etal).^(-1);
%     
%     % Get stress -------------------------------------------------------- %
%     D.tauss     = 2.*D.etass.*E;                        % [ Pa ]
    
    % Grain size -------------------------------------------------------- %
    Rss         = G.Kam97{4}.*(D.tauII.*ssc).^(-G.Kam97{5});  % [ m ]
    
    Rss         = -0.5.*...
        (-G.Kam97{3}.*D.epstot.*esc.*Rss - ...
        sqrt(G.Kam97{3}^2.*(D.epstot.*esc).^2.*Rss.^2 + ...
        4*G.Kam97{3}.*D.epstot*esc.*k))./(G.Kam97{3}.*D.epstot.*esc);
    
    D.d         =   D.d.*G.itfac + Rss.*(1-G.itfac);
%     D.Rinf      = D.Rinf.*0.9 + Rss.*(1-0.9);
    
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
%             semilogy(i,D.Rinf(:),'k.',i,Rss,'ro')
%             title('Grain Size')
%             hold on
%             pause(.1)
%     end
% end
% ======================================================================= %
% ============================ END ====================================== %
% ======================================================================= %
end
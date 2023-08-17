function [D] = GSE_SS_BRAUN99(D,G,GC,T)
% ======================================================================= %
%                                                                         %
% ======================================================================= %

D.lab       = 'Br99a';

% Grain growth ---------------------------------------------------------- %
switch lower(G.growth)
    case 'yes'
        D.lab   = 'Br99b';
        k       = G.Br99{1}.*exp(-G.Br99{2}./GC.RG./T.T);      % [ mm^2/s ]
end
% ----------------------------------------------------------------------- %

% Scale for strain rate and stress for second invariant use ------------- %
esc         = 2/sqrt(3);
ssc         = sqrt(3);
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% Assuming both mechanisms act simultaneously
% for i=1:200
%     % Dislocation creep viscosity --------------------------------------- %
%     D.etal = 0.5*R.Al^(-1/R.n)*...
%         exp(R.Ql/GC.RG./T./R.n).*D.epsss.^(1/R.n-1);    % [ Pa s ]
%
%     % Diffusion creep viscosity ----------------------------------------- %
%     % Grain size needs to be in m! -------------------------------------- %
%     D.etaf = 0.5*R.Af^(-1)*...
%         exp(R.Qf/GC.RG./T).*D.Rinf.^R.m;                % [ Pa s ]
%
%     % Effective viscosity ----------------------------------------------- %
%     D.etass     = (1./D.etaf + 1./D.etal).^(-1);
%
%     % Get stress -------------------------------------------------------- %
%     D.tauss     = 2.*D.etass.*E;                        % [ Pa ]

% Grain size -------------------------------------------------------- %
switch lower(G.growth)
    case 'no'
        Rinf        = G.Br99{4}.*(D.tauII.*ssc).^(-G.Br99{5});  % [ m ]
    case 'yes'
        gammaG      = D.etal./(D.etal + D.etaf);
        Rinf        = G.Br99{4}.*(D.tauII.*ssc).^(-G.Br99{5});  % [ m ]
        Rinf        = 0.5.*(D.epstot.*esc.*Rinf + ...
            sqrt((D.epstot.*esc).^2.*Rinf.^2 + ...
            4.*k.*D.epstot.*esc.*gammaG.*G.Br99{3}))./(D.epstot.*esc);
end

D.d     =   D.d.*G.itfac + Rinf.*(1-G.itfac);

%     D.efss      = D.tauss./2./D.etaf;
%     D.epsss     = D.tauss./2./D.etal;
%     D.etotss    = D.efss + D.epsss;

% switch N.debug
%     case 3
%         if (max(size(D.Rinf))>1)
%             error('ERROR! Dimension for debug check to large to plot!')
%         end
%         figure(1)
%         subplot(2,2,1)
%         plot(i,log10(D.etass(:)),'k+',...
%             i,log10(etaf(:)),'b.',...
%             i,log10(etal(:)),'ro')
%         title('Viscosity')
%         hold on
%         subplot(2,2,2)
%         plot(i,log10(D.tauss(:)),'k.')
%         title('Stress')
%         hold on
%         subplot(2,2,3)
%         plot(i,log10(D.efss(:)),'b.',...
%             i,log10(D.epsss(:)),'ro',...
%             i,log10(D.etotss(:)),'k+')
%         title('Strain Rate')
%         hold on
%         subplot(2,2,4)
%         semilogy(i,D.Rinf(:),'k.',i,Rinf,'ro')
%         title('Grain Size')
%         hold on
%         pause(.1)
% end
% end
% ======================================================================= %
% ============================ END ====================================== %
% ======================================================================= %
end
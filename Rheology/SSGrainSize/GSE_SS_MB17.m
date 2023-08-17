function [D] = GSE_SS_MB17(D,G,GC,R,T)
% All constants are used in SI units
%                                                                         %
% ======================================================================= %

D.lab       = 'MB17';

D.d      = D.d.*sqrt(pi/2);        % Initial interface roughness

% Interface roughness growth -------------------------------------------- %
GI      = G.MB17{1}.*exp(-G.MB17{2}./(GC.RG.*T.T)).*...
    G.MB17{5};                  % [ m^p s^-1 ]
% ----------------------------------------------------------------------- %

% Rheological compliances ----------------------------------------------- %
A           = R.Al.*exp(-R.El./GC.RG./T.T);
B           = R.Af.*exp(-R.Ef./GC.RG./T.T);
% ----------------------------------------------------------------------- %

% Calculate work partitioning ------------------------------------------- %
dE      = R.El-G.MB17{2};               % [ J/mol ]
F       = 3.*G.MB17{9}.*G.MB17{7}.*G.MB17{1}./...
    2./G.MB17{10}./R.Al;                  % [ Pa^(n+1) m^(p+1) ]
Cf      = 35 - 9e-5.*dE - 0.5.*log10(F);
kf      = 1.7e-19.*dE.^(3.4).*F.^(0.081);

fIp     = G.MB17{11}.*exp(-Cf.*(T.T./1000).^kf);
% ----------------------------------------------------------------------- %

% Constants ------------------------------------------------------------- %
CI      = G.MB17{6}.*GI;                % [ m^p s^-1 ]
% ----------------------------------------------------------------------- %

% D.epsss     = E;

% % Assuming both mechanisms act simultaneously---------------------------- %
% for i=1:100
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
    
    psi         = 2.*D.tauII.*D.epsln;
    
    % Grain size -------------------------------------------------------- %
    Rf          = (B./A./D.tauss.^(R.n-1)).^(1/R.m);
    
    Rm          = G.MB17{8}.*Rf;
    
    fI          = fIp.*...
        (Rm.^G.MB17{12}./(Rm.^G.MB17{12} + D.rinf.^G.MB17{12}));
    
    DI          = fI./G.MB17{7}/G.MB17{6};
    
    rinf        = (CI./G.MB17{4}./DI).^(1/(G.MB17{4}+1)).*...
        (1./psi).^(1/(G.MB17{4}+1));
    
    D.rinf      = D.rinf.*0.9 + rinf.*(1-0.9);
    
    D.Rinf      = D.rinf./sqrt(pi/2);
    
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
function [D] = EffectiveViscosity(D,R)
% ======================================================================= %
%                                                                         %
% ======================================================================= %

switch lower(R.cetaeff)
    case 'harmonic'
        D.eta   =   (1./D.etal + 1./D.etaf).^(-1);
    case 'minimum'
        D.eta   =   min(D.etal,D.etaf);
    otherwise
        D.eta   =   D.etaf;
end

switch lower(R.plast)
    case {'drucker-prager','yes'}
        D.eta   =   min(D.eta,D.etay);
end

end
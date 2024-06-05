function [D] = SteadyStateGrainSize(D,R,G,GC,T)
% ======================================================================= %
%                                                                         %
% ======================================================================= %

addpath('SSGrainSize/')

switch lower(G.GSE)
    case {'br99a','br99b'}
        D       =   GSE_SS_BRAUN99(D,G,GC,T);
    case {'ae07','be09'}
        D       =   GSE_SS_AE07(D,G,GC,R,T);
    case 'dan17'
        D       =   GSE_SS_DAN17(D,G,GC,R,T);
    case 'kam97'
        D       =   GSE_SS_KAM97(D,GC,G,T);
    case 'ge20'
        D       =   GSE_SS_Ge20(D,G,GC,T);
    case 'ro11'
        D       =   GSE_SS_ROZ11(D,G,GC,R,T);
end

rmpath('SSGrainSize/')

end
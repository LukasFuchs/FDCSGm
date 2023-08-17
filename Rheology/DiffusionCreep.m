function D = DiffusionCreep(D,R,T,GC)

switch lower(R.pressure)
    case 'yes'
        R.Qf    =   R.Ef + D.P.*R.Vf;
end

D.etaf      =   0.5.*R.Af.^(-1).*...
    exp(R.Qf./GC.RG./T.T).*D.d.^R.m;     % [ Pa s ]

end
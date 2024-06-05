function D = DislocationCreep(D,R,T,GC)

switch lower(R.pressure)
    case 'yes'
        R.Ql    =   R.El + D.P.*R.Vl;
end

D.etal  =   0.5.*R.Al.^(-1./R.n).*...
    exp(R.Ql./GC.RG./T.T./R.n).*D.epsl.^(1./R.n-1);    % [ Pa s ]

end
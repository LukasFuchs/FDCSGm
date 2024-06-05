function T = Diffusion(B,T,Q,rho,Py,dt,N)
switch Py.scale
    case 'no'
        switch lower(B.DiffMethod)
            case 'explicit'
                T   =   SolveDiff2Dexplicit(T,Q,rho,dt,Py,N,B);
            case 'implicit'
                T   =   SolveDiff2Dimplicit(T,Q,rho,dt,Py,N,B);
            case 'adi'
                T   =   SolveDiff2DADI(T,Q,rho,dt,Py,N,B);
            case 'cnv'
                T   =   SolveDiff2DCNV(T,Q,rho,dt,Py,N,B);
            case 'none'
            otherwise
                error('Diffusion scheme is not defined!')
        end
    case 'yes'
        switch lower(B.DiffMethod)
            case 'explicit'
                T     =   SolveDiff2DexplicitSc(T,Q,dt,N,B);
            case 'implicit'
                T     =   SolveDiff2DimplicitSc(T,Q,dt,N,B);
            case 'adi'
                T     =   SolveDiff2DADISc(T,Q,dt,N,B);
            case 'none'
            otherwise
                error('Diffusion scheme is not defined!')
        end
end
end
function [answer,T] = CheckBreakCriteria(it,T,D,M,Pl,ID,Py)

answer      = 'no';

switch lower(Py.scale)
    case 'yes'
        T.tind  =   find(T.time(1:it) >= (T.time(it)-0.02),1,'first');
        epsC    =   1e-2;
    case 'no'
        tspan   =   50;  % [ Ma ]
        tspan   =   tspan*1e6*(60*60*24*365.25);    % [ s ]
        T.tind  =   find(T.time(1:it) >= (T.time(it)-tspan),1,'first');
        epsC    =   1e-16; 
end

if (T.tind > 1)
    epsV    =   std(D.meanV(T.tind:it));
else    
    epsV    =   1e3;
end

if (mod(it,10)==0||it==1)
    disp(['eps_V = ',sprintf('%g',epsV)])
end

if (T.time(it) >= T.tmax)
    if strcmp(getenv('OS'),'Windows_NT')
        set(figure(1),'position',[1.8,1.8,766.4,780.8]);
    else
        set(figure(1),'position',[1.8,26,866.2,949]);
    end
    disp('Maximale Zeit erreicht. Zeitschleife unterbrochen')
    T.indtime   = find(T.time(2:end)==0,1,'first');
    figure(1) % ----------------------------------------------------- %
    clf
    switch Py.eparam
        case 'const'
            ax1=subplot(2,1,1);
            plotfield(D.T,M.X,M.Z,Pl,'contourf',...
                '$$T$$','quiver',ID.vx,ID.vz)
            colormap(ax1,flipud(Pl.lajolla))
            ax2=subplot(2,1,2);
            plotfield(ID.v,M.X,M.Z,Pl,'pcolor',...
                '$$v$$')
            colormap(ax2,Pl.imola);
        case 'variable'
            ax1=subplot(2,1,1);
            plotfield(D.T,M.X,M.Z,Pl,'contourf',...
                '$$T$$','quiver',ID.vx,ID.vz)
            colormap(ax1,flipud(Pl.lajolla))
            ax2=subplot(2,1,2);
            plotfield(log10(D.eta),M.X,M.Z,Pl,'contourf',...
                '$$\eta$$','quiver',ID.vx,ID.vz)
            colormap(ax2,flipud(Pl.lapaz))
    end
    answer = 'yes';
elseif( epsV <= epsC && it > 50 )
    if strcmp(getenv('OS'),'Windows_NT')
        set(figure(1),'position',[1.8,1.8,766.4,780.8]);
    else
        set(figure(1),'position',[1.8,26,866.2,949]);
    end
%     h           =   figure(1);
    T.indtime   = find(T.time(2:end)==0,1,'first');
    disp('Konvektion erreicht steady state.')
    figure(1) % ----------------------------------------------------- %
    clf
    switch Py.eparam
        case 'const'
            ax1=subplot(2,1,1);
            plotfield(D.T,M.X,M.Z,Pl,'contourf',...
                '$$T$$','quiver',ID.vx,ID.vz)
            colormap(ax1,flipud(Pl.lajolla))
            ax2=subplot(2,1,2);
            plotfield(ID.v,M.X,M.Z,Pl,'pcolor',...
                '$$v$$')
            colormap(ax2,Pl.imola)
        case 'variable'
            ax1=subplot(2,1,1);
            plotfield(D.T,M.X,M.Z,Pl,'contourf',...
                '$$T$$','quiver',ID.vx,ID.vz)
            colormap(ax1,flipud(Pl.lajolla))
            ax2=subplot(2,1,2);
            plotfield(log10(D.eta),M.X,M.Z,Pl,'contourf',...
                '$$\eta$$','quiver',ID.vx,ID.vz)
            colormap(ax2,flipud(Pl.lapaz))
    end
    answer = 'yes';
elseif(it == T.itmax)
    T.indtime   = T.itmax;
    disp('Maximale Anzahl der Iterationen erreicht.')
end

end

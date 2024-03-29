function Pl = PlotData_dim(it,Pl,T,D,M,ID,Py)
Pl.time     =   ...
    [{['@ Iteration: ',sprintf('%i',it)]};...
    {['Time: ',sprintf('%2.2e',T.time(it)/(1e6*(365.25*24*60*60))),' [Ma]']}];
if (mod(it,Pl.tstpinc)==0||it==1)
    switch Pl.plotfields
        case 'yes'
            figure(1)
            clf
            switch Py.eparam
                case 'const'
                    ax1=subplot(2,1,1);
                    plotfield(D.T,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                        'T [K]','quiver',ID.vx,ID.vz)
                    caxis([Py.T0 Py.T0+Py.DeltaT])
                    colormap(ax1,Pl.lajolla)
                    ax2=subplot(2,1,2);
                    plotfield(ID.v,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
                        'v [m/s]')
                    colormap(ax2,Pl.imola)
                case 'variable'
                    ax1=subplot(2,1,1);
                    plotfield(D.T,M.X./1e3,M.Z./1e3,Pl,'contourf',...
                        'T [K]','quiver',ID.vx,ID.vz)
                    colormap(ax1,Pl.lajolla)
                    ax2=subplot(2,1,2);
                    plotfield(log10(D.eta),M.X./1e3,M.Z./1e3,Pl,'contourf',...
                        '$$\eta$$ [Pa s]','quiver',ID.vx,ID.vz)
                    colormap(ax2,flipud(Pl.lapaz))
            end
            switch Pl.savefig
                case 'yes'
                    saveas(figure(1),...
                        [M.ModDir,'/Field',num2str(it)],'png')
                    % Capture the plot as an image
                    frame       = getframe(Pl.h);
                    im          = frame2im(frame);
                    [imind,cm]  = rgb2ind(im,256);
                    % Write to the GIF File
                    if it == 1
                        imwrite(imind,cm,Pl.filename,'gif', 'Loopcount',inf);
                    else
                        imwrite(imind,cm,Pl.filename,'gif','WriteMode','append');
                    end
            end
        case 'no'
            disp(['Iteration: ',sprintf('%i',it)])
    end
end
end
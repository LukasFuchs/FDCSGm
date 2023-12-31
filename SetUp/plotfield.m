function plotfield(field,X,Z,Pl,plparam,titparam,plparam2,qx,qz)

if nargin < 7
    plparam2    =   'none';
end

fmin    = min(min(field));
fmax    = max(max(field));

if fmin == fmax
    fmax = fmin + 1;
end

switch lower(plparam)
    case 'pcolor'
        pcolor(X,Z,field)
    case 'contour'
        contour(X,Z,field,10)
    case 'contourf'
        contourf(X,Z,field,10)
    case 'scatter'
        scatter(X,Z,4,field,'filled')
    otherwise
        error('Plot Method unknown! Check plotfield call')
end

switch lower(plparam2)
    case 'quiver'
        hold on
        quiver(X(1:Pl.inc:end,1:Pl.inc:end),...
            Z(1:Pl.inc:end,1:Pl.inc:end),...
            qx(1:Pl.inc:end,1:Pl.inc:end),...
            qz(1:Pl.inc:end,1:Pl.inc:end),1,'w','LineWidth',2)
        hold off
    case 'contour'
        hold on
        contour(X,Z,qx,10,'k','LineWidth',2,'ShowText','on')
        hold off
    case 'contour2'
        hold on
        contour(X,Z,qx,[qz qz],'k','LineWidth',1)
        hold off
    case 'contoury'
        hold on
        contour(X,Z,qx,10,'y--')
        hold off
    case 'scatter'
        hold on
        scatter(X,Z,1,qx,'filled')
        hold off
end

xlabel(Pl.xlab,'Interpreter','latex')
ylabel(Pl.zlab,'Interpreter','latex')
if ~isfield(Pl,'time')
    title(titparam,'Interpreter','latex')
else
    if ~isempty(titparam)
        title([titparam;Pl.time],'Interpreter','latex')
    else
        title(Pl.time,'Interpreter','latex')
    end
end
cb = colorbar; cb.TickLabelInterpreter = 'latex';
if isfield(Pl,'cbtitle')
    title(cb,Pl.cbtitle,'Interpreter','latex')
end
shading interp, lighting phong
axis normal; axis equal, axis tight; %caxis([fmin fmax])
axis([min(X(:)) max(X(:)) min(Z(:)) max(Z(:))])
set(gca,'FontWeight','Bold','FontSize',12,...
    'TickLabelInterpreter','latex')
box on

end
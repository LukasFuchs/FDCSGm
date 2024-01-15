function PlotProfile(X,Y,tit,xlab,ylab,xmin,xmax,ymin,ymax)

if nargin < 3
    param1 = 0; 
else
    param1 = 1; 
end

if nargin < 6
    param2 = 0;
else
    param2 = 1; 
end

plot(X,Y,'LineWidth',2)
set(gca,'FontWeight','Bold','FontSize',12,'LineWidth',2,...
    'TickLabelInterpreter','latex')
box on; 
if param1 == 1
   title(tit,'Interpreter','latex') 
   xlabel(xlab,'Interpreter','latex')
   ylabel(ylab,'Interpreter','latex')
end
if param2 == 1
    axis([xmin xmax ymin ymax])
end
    

end


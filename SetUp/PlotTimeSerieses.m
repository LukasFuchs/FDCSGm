function PlotTimeSerieses(Py,T,D,M,N,Ger)

if nargin < 6
    Ger = 0;
    Benchmark = 0;
    f = 3;
else
    Benchmark = 1;
    f = 4;
end

if Benchmark
    switch Py.eparam
        case 'const'
            if (Py.eta0 == 1e23)
                blmod   = 1;
            elseif (Py.eta0 == 1e22)
                blmod   = 2;
            elseif (Py.eta0 == 1e21)
                blmod   = 3;
            end
        case 'variable'
            if Py.c == 0
                blmod   = 4;
            else
                blmod   = 5;
            end
    end
end

set(figure(2),'position',[1,1,1536,788.8]);

figure(2)
clf
subplot(f,1,1)
plot(T.time(1:T.indtime),D.Nus(1:T.indtime),...
    'LineWidth',2)
if Benchmark
    hold on
    plot(T.time(1:T.indtime),ones(T.indtime,1).*Ger(1,blmod),'r--')
end
set(gca,'FontWeight','Bold');
xlabel('t'); ylabel('Nus')
title('Nusselt Number')

subplot(f,1,2)
plot(T.time(1:T.indtime),D.meanV(1:T.indtime),...
    'LineWidth',2)
if Benchmark
    hold on
    %     plot(T.time(1:T.indtime),D.meanV2(1:T.indtime),'r--',...
    %     'LineWidth',2)
    plot(T.time(1:T.indtime),ones(T.indtime,1).*Ger(2,blmod),'r--')
end
set(gca,'FontWeight','Bold');
xlabel('t'); ylabel('VRMS')
title('Root Mean Square Velocity')


subplot(f,1,3)
if Benchmark
    plot(D.T(:,round(N.nx1/2)),M.z,'LineWidth',2)
    hold on
    plot(Ger(3,blmod),M.H+Ger(4,blmod),'sk')
    plot(Ger(5,blmod),M.H+Ger(6,blmod),'sk')
else
    plot(D.meanT(:,T.indtime),M.z,'LineWidth',2)
end
set(gca,'FontWeight','Bold');
xlabel('T'); ylabel('z')
title('Temperature Profile')

if Benchmark
    subplot(4,1,4)
    plot([D.dTtop(1),D.dTtop(end),D.dTbot(end),D.dTbot(1)],'o',...
        'MarkerFaceColor','k')
    if Benchmark
        hold on
        plot([Ger(7,blmod),Ger(8,blmod),Ger(9,blmod),Ger(10,blmod)],'s',...
            'MarkerFaceColor','r')
    end
    set(gca,'FontWeight','Bold','yscale','log');
    xlabel('\DeltaT'); ylabel('')
    title('Temperature difference at corners')
    legend('Model','Benchmark')
end

end

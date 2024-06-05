% ----------------------------------------------------------------------- %
% Visualize data for a viscous inclusion in an incompressible medium for  %
% different orientations and viscosity contrasts                          %
%                                                                         %
% ----------------------------------------------------------------------- %
% LF - 25.07.2023 -                                                       %
% ======================================================================= %
clear
clc
clf

DDir        =   'data/Ell_a_300_b_100/';

DefType     =   'PureShear';

angle       =   [0 22.5 45 90];

DeltaE      = zeros(40,length(angle));
psi         = zeros(40,length(angle));
eII         = zeros(40,length(angle));
tauII       = zeros(40,length(angle));
psima       = zeros(40,length(angle));
eIIma       = zeros(40,length(angle));
tauIIma     = zeros(40,length(angle));
psimi       = zeros(40,length(angle));
eIImi       = zeros(40,length(angle));
tauIImi     = zeros(40,length(angle));
psistd      = zeros(40,length(angle));
eIIstd      = zeros(40,length(angle));
tauIIstd    = zeros(40,length(angle));

psimat      = zeros(40,length(angle));
eIImat      = zeros(40,length(angle));
tauIImat    = zeros(40,length(angle));
psimatma    = zeros(40,length(angle));
eIImatma    = zeros(40,length(angle));
tauIImatma  = zeros(40,length(angle));
psimatmi    = zeros(40,length(angle));
eIImatmi    = zeros(40,length(angle));
tauIImatmi  = zeros(40,length(angle));
psimatstd   = zeros(40,length(angle));
eIImatstd   = zeros(40,length(angle));
tauIImatstd = zeros(40,length(angle));

col         = colormap(lines(length(angle)));
legendinfo  = cell(length(angle),1);

set(figure(1),'Position',[1.8,1.8,766.4,780.8]);
set(figure(2),'Position',[769.8,1.8,766.4,780.8]);
% set(figure(3),'Position',[1.8,1.8,766.4,780.8]);
% set(figure(4),'Position',[1.8,1.8,766.4,780.8]);

for k = 1:length(angle)
    % data = [eta2'./Py.eta0 psiinc1 eIIinc tauIIinc];
    DATA    = load([DDir,'/Data_',DefType,'_',num2str(angle(k)),'.mat']);

    DeltaE(:,k)     =   DATA.data(:,1);
    psi(:,k)        =   DATA.data(:,2);
    psima(:,k)      =   DATA.data(:,3);
    psimi(:,k)      =   DATA.data(:,4);
    psistd(:,k)     =   DATA.data(:,5);

    eII(:,k)        =   DATA.data(:,6);
    eIIma(:,k)      =   DATA.data(:,7);
    eIImi(:,k)      =   DATA.data(:,8);
    eIIstd(:,k)     =   DATA.data(:,9);

    tauII(:,k)      =   DATA.data(:,10);
    tauIIma(:,k)    =   DATA.data(:,11);
    tauIImi(:,k)    =   DATA.data(:,12);
    tauIIstd(:,k)   =   DATA.data(:,13);

    psimat(:,k)         =   DATA.data(:,14);
    psimatma(:,k)       =   DATA.data(:,15);
    psimatmi(:,k)       =   DATA.data(:,16);
    psimatstd(:,k)      =   DATA.data(:,17);

    eIImat(:,k)         =   DATA.data(:,18);
    eIImatma(:,k)       =   DATA.data(:,19);
    eIImatmi(:,k)       =   DATA.data(:,20);
    eIImatstd(:,k)      =   DATA.data(:,21);

    tauIImat(:,k)       =   DATA.data(:,22);
    tauIImatma(:,k)     =   DATA.data(:,23);
    tauIImatmi(:,k)     =   DATA.data(:,24);
    tauIImatstd(:,k)    =   DATA.data(:,25);

    clear DATA

    figure(1)
    subplot(3,2,1)
    hold on       
    errorbar(DeltaE(:,k),eII(:,k),eIIstd(:,k),'LineWidth',2,'Marker','s')    
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'xscale','log','yscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')
    title([DefType,' Inclusion'])
    subplot(3,2,2)
    hold on    
    errorbar(DeltaE(:,k),eIImat(:,k),eIImatstd(:,k),'LineWidth',2,'Marker','s')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'xscale','log','yscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')   
    title([DefType,' Matrix'])

    subplot(3,2,3)
    hold on
    errorbar(DeltaE(:,k),tauII(:,k),tauIIstd(:,k),'LineWidth',2,'Marker','s') 
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'xscale','log','yscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \tau_{II} \rangle [s^{-1}]$$','Interpreter','latex')    
    subplot(3,2,4)
    hold on
    errorbar(DeltaE(:,k),tauIImat(:,k),tauIImatstd(:,k),'LineWidth',2,'Marker','s')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'xscale','log','yscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \tau_{II} \rangle [s^{-1}]$$','Interpreter','latex')

    subplot(3,2,5)
    hold on
    errorbar(DeltaE(:,k),psi(:,k),psistd(:,k),'LineWidth',2,'Marker','s')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'xscale','log','yscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \psi_{II} \rangle [s^{-1}]$$','Interpreter','latex')
    subplot(3,2,6)
    hold on
    errorbar(DeltaE(:,k),psimat(:,k),psimatstd(:,k),'LineWidth',2,'Marker','s')
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'xscale','log','yscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
    ylabel('$$\langle \psi_{II} \rangle [s^{-1}]$$','Interpreter','latex')

    figure(2)    
%     subplot(2,1,1)
    hold on
    p(k) = plot(eII(:,k),tauII(:,k),'LineWidth',2,'Color',col(k,:));
    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
        'yscale','lin','xscale','lin','TickLabelInterpreter','latex');
    box on
    xlabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')
    ylabel('$$\langle \tau_{II} \rangle [MPa]$$','Interpreter','latex')
    legendinfo{k} = ['\theta = ',num2str(angle(k))];
    %     axis([1e-16 1e-10 1e-2 1e3])
    title([DefType,' Inclusion'])

%     subplot(2,1,2)
%     hold on
%     plot(eIImat(:,k),tauIImat(:,k),'LineWidth',2,'Color',col(k,:));    
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'yscale','log','xscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     ylabel('$$\langle \tau_{II} \rangle [MPa]$$','Interpreter','latex')
%     title([DefType,' Matrix'])

%     figure(3)
%     subplot(3,2,1)
%     hold on    
%     plot(DeltaE(:,k),eIIma(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     title([DefType,' Inclusion Max'])    
%     subplot(3,2,2)
%     hold on    
%     plot(DeltaE(:,k),eIImatma(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')   
%     title([DefType,' Matrix Max'])
% 
%     subplot(3,2,3)
%     hold on
%     plot(DeltaE(:,k),tauIIma(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \tau_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     subplot(3,2,4)
%     hold on
%     plot(DeltaE(:,k),tauIImatma(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \tau_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     
%     subplot(3,2,5)
%     hold on
%     plot(DeltaE(:,k),psima(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \psi_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     subplot(3,2,6)
%     hold on
%     plot(DeltaE(:,k),psimatma(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \psi_{II} \rangle [s^{-1}]$$','Interpreter','latex')
% 
%     figure(4)
%     subplot(3,2,1)
%     hold on    
%     plot(DeltaE(:,k),eIImi(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     title([DefType,' Inclusion Min'])    
%     subplot(3,2,2)
%     hold on    
%     plot(DeltaE(:,k),eIImatmi(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \epsilon_{II} \rangle [s^{-1}]$$','Interpreter','latex')   
%     title([DefType,' Matrix Min'])
% 
%     subplot(3,2,3)
%     hold on
%     plot(DeltaE(:,k),tauIImi(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \tau_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     subplot(3,2,4)
%     hold on
%     plot(DeltaE(:,k),tauIImatmi(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \tau_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     
%     subplot(3,2,5)
%     hold on
%     plot(DeltaE(:,k),psimi(:,k),'LineWidth',2,'Color',col(k,:))
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \psi_{II} \rangle [s^{-1}]$$','Interpreter','latex')
%     subplot(3,2,6)
%     hold on
%     plot(DeltaE(:,k),psimatmi(:,k),'LineWidth',2,'Color',col(k,:))    
%     set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',13,...
%         'xscale','log','yscale','log','TickLabelInterpreter','latex');
%     box on
%     xlabel('$$\Delta ( \eta )$$','Interpreter','latex');
%     ylabel('$$\langle \psi_{II} \rangle [s^{-1}]$$','Interpreter','latex')
end

legend(p,legendinfo,'Location','SouthWest')


function ID = CheckContinuum(ID,N,M,Ma,Pl)

if strcmp(getenv('OS'),'Windows_NT')
    set(figure(76),'position',[1.0,1.0,1536.0,788.8]);
else
    set(figure(76),'position',[-1919,1,960,988]);
end

if ~isempty(Ma)
    XM = reshape(Ma.XM,[N.nmz*N.nz1,N.nmx*N.nx1]);
    ZM = reshape(Ma.ZM,[N.nmz*N.nz1,N.nmx*N.nx1]);
    CM = reshape(Ma.C,[N.nmz*N.nz1,N.nmx*N.nx1]);
end

Pl.cbtitle  =   {};
ID      =   VelGradTensor(ID,N);
figure(76)
clf
ax1 = subplot(2,3,1);
plotfield(ID.dvxdx,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '$$\frac{\partial v_x}{\partial x}$$')
caxis([-1e-14 1e-14])
colormap(ax1,Pl.vik)
ax2 = subplot(2,3,2);
plotfield(ID.dvxdz,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '$$\frac{\partial v_x}{\partial z}$$')
caxis([-1e-14 1e-14])
colormap(ax2,Pl.vik)
ax3 = subplot(2,3,4);
plotfield(ID.dvzdx,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '$$\frac{\partial v_z}{\partial x}$$')
caxis([-1e-14 1e-14])
colormap(ax3,Pl.vik)
ax4 = subplot(2,3,5);
plotfield(ID.dvzdz,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '$$\frac{\partial v_z}{\partial z}$$')
caxis([-1e-14 1e-14])
colormap(ax4,Pl.vik)
ax5 = subplot(2,3,3);
plotfield(ID.div,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '$$div(v)$$')
caxis([-1e-20 1e-20])
if ~isempty(Ma)
    hold on
    ind1    =   CM==1;
    scatter1 = scatter(XM(ind1)./1e3,ZM(ind1)./1e3,1,'MarkerFaceColor','k');
    alpha(scatter1,0.1)
end
colormap(ax5,Pl.vik)
ax6 = subplot(2,3,6);
plotfield(ID.rot,M.X./1e3,M.Z./1e3,Pl,'pcolor',...
    '$$rot(v)$$')
caxis([-1e-14 1e-14])
if ~isempty(Ma)
    hold on
    ind1    =   CM==1;
    scatter1 = scatter(XM(ind1)./1e3,ZM(ind1)./1e3,1,'MarkerFaceColor','k');
    alpha(scatter1,0.1)
end
colormap(ax6,Pl.vik)


end
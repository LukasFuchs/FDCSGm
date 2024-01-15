function [L1PrErr2,L1VxErr2,dx] = M2Di_Linear(resol,solver, SuiteSparse,noisy)
% M2Di Linear Stokes: Thibault Duretz, Ludovic Raess, Yuri Podladchikov  - Unil 2017
% Slightly modified to fit to the problem. 
% ----------------------------------------------------------------------- %
%% Physics
Lx      = 6;                                                                % Box width
Ly      = 6;                                                                % Box height
mus0    = 1;                                                                % Background viscosity
mus_i   = 1e4;                                                              % Inclusion viscosity
rad     = 1;                                                                % Inclusion radius
%% Numerics
nx      = resol;                                                            % Grid points in x
ny      = nx;                                                               % Grid points in y
gamma   = 1e3;                                                              % Numerical compressibility
tol_lin = 1e-12;                                                            % Linear solver tolerance
niter   = 10;                                                               % Max number of linear iterations
%% Preprocessing
tic
dx = Lx/nx;                                                                 % cell size in x
dy = Ly/ny;                                                                 % cell size in y
xv = -Lx/2:dx:Lx/2;           yv = -Ly/2:dy:Ly/2;           [XV2, YV2] = ndgrid(xv, yv); % cell coord. grid
xc = -Lx/2+dx/2:dx:Lx/2-dx/2; yc = -Ly/2+dy/2:dy:Ly/2-dy/2; [XC2, YC2] = ndgrid(xc, yc); % cell coord. grid
%% Initial conditions
Pt                                =      zeros(nx  ,ny  );                  % Pressure
Vx                                =      zeros(nx+1,ny  );                  % x velocity
Vy                                =      zeros(nx  ,ny+1);                  % y velocity
mus                               =  mus0*ones(nx  ,ny  );                  % viscosity at cell center
musv                              =  mus0*ones(nx+1,ny+1);                  % viscosity at vertex
mus( sqrt(XC2.^2 + YC2.^2) < rad) = mus_i;                                  % define inclusion
musv(sqrt(XV2.^2 + YV2.^2) < rad) = mus_i;
mus_axy                           = musv(2:end-1, 2:end-1);
%% Numbering Pt and Vx,Vy
NumVx  = reshape(1:(nx+1)*ny,nx+1,ny  );
NumVy  = reshape(1:nx*(ny+1),nx  ,ny+1); NumVyG = NumVy + max(NumVx(:));    % G stands for Gobal numbering
NumPt  = reshape(1:nx*ny    ,nx  ,ny  ); NumPtG = NumPt + max(NumVyG(:));   % G stands for Gobal numbering
cpu(1)=toc;
%% DANI analytical solution
tic
[ Vx_W, Vx_E, Vx_S, Vx_N, Vy_W, Vy_E, Vy_S, Vy_N, Pa, Vxa, Vya ] = Dani_Solution_vec(xv,yv,xc,yc,rad,mus_i ,nx,ny);
cpu(2)=toc; display(['Time Dani sol     = ', num2str(cpu(2))]);
%% Boundary Conditions on velocities [W E S N]
tic
ibcVxW = NumVx(1,:)';         ibcVxE = NumVx(end,:)';                       % Indexes of BC nodes for Vx East/West
ibcVxS = NumVx(2:end-1,1);    ibcVxN = NumVx(2:end-1,end);                  % Indexes of BC nodes for Vx North/South
ibcVyS = NumVyG(:,1);         ibcVyN = NumVyG(:,end);                       % Indexes of BC nodes for Vy North/South
ibcVyW = NumVyG(1,2:end-1)';  ibcVyE = NumVyG(end,2:end-1)';                % Indexes of BC nodes for Vy East/West
ibc    = [ ibcVxW; ibcVxE; ibcVyS; ibcVyN ];                                % Group all indexes
ibcNC  = [ ibcVxS; ibcVxN; ibcVyW; ibcVyE ];                                % Non Confroming to the physical boundary
vBc    = [ Vx_W;          Vx_E;          Vy_S;          Vy_N         ];     % Group all values
vBcNC  = [ Vx_S(2:end-1); Vx_N(2:end-1); Vy_W(2:end-1); Vy_E(2:end-1)];     % Non Confroming to the physical boundary
cpu(3)=toc;
%% grad and div blocs
tic
iVxC   = NumVx;                                                             % dP/dx
iPtW   =  ones(size(iVxC));    iPtW(1:end-1,:) = NumPt;
iPtE   =  ones(size(iVxC));    iPtE(2:end  ,:) = NumPt;
cPtW   =  ones(size(iVxC))/dx; cPtW([1 end],:) = 0;
cPtE   = -ones(size(iVxC))/dx; cPtE([1 end],:) = 0;
Idx    = [ iVxC(:); iVxC(:) ]';
Jdx    = [ iPtW(:); iPtE(:) ]';
Vdx    = [ cPtW(:); cPtE(:) ]';
iVyC   = NumVyG;                                                            % dP/dy
iPtS   =  ones(size(iVyC));    iPtS(:,1:end-1) = NumPt;
iPtN   =  ones(size(iVyC));    iPtN(:,2:end  ) = NumPt;
cPtS   =  ones(size(iVyC))/dy; cPtS(:,[1 end]) = 0;
cPtN   = -ones(size(iVyC))/dy; cPtN(:,[1 end]) = 0;
Idy    = [ iVyC(:); iVyC(:) ]';
Jdy    = [ iPtS(:); iPtN(:) ]';
Vdy    = [ cPtS(:); cPtN(:) ]';
%% PP block
iPt    = NumPt;  I = iPt(:)';  J = I;                                       % Eq. index center (pressure diagonal)
V      = ones(nx*ny,1)./gamma;                                              % Center coeff.
if SuiteSparse==1, PP = sparse2(I,J,V); else                                % Matrix assembly
                   PP =  sparse(I,J,V); end
%% Block UU
mus_W  = zeros(size(Vx));  mus_W(2:end  , :     ) = mus;                    % Viscosities (W,E,S,N)
mus_E  = zeros(size(Vx));  mus_E(1:end-1, :     ) = mus;
mus_S  = zeros(size(Vx));  mus_S(2:end-1,1:end-1) = mus_axy;
mus_N  = zeros(size(Vx));  mus_N(2:end-1,2:end  ) = mus_axy;
iVx    = NumVx;                                                             % Eq. index for Vx (C,W,E,S,N)
iVxW   = NumVx(1:end-1, :     );
iVxE   = NumVx(2:end  , :     );
iVxS   = NumVx(2:end-1,1:end-1);
iVxN   = NumVx(2:end-1,2:end  );
cVxC   = (mus_W+mus_E)/dx/dx + (mus_S+mus_N)/dy/dy;                         % Center coeff.
scVx   = max(cVxC(:));                   cVxC([1,end],:) = scVx;            % Scaling factor for Vx Dirichlet values
cVxW   = -mus_W(2:end  , :     )/dx/dx;  cVxW([1,end],:) = 0;               % West coeff.
cVxS   = -mus_S(2:end-1,1:end-1)/dy/dy;                                     % South coeff.
Iuu    = [  iVx(:); iVxE(:); iVxN(:) ]';                                    % Triplets [I,J,V]
Juu    = [  iVx(:); iVxW(:); iVxS(:) ]';
Vuu    = [ cVxC(:); cVxW(:); cVxS(:) ]';
%% Block VV
mus_W  = zeros(size(Vy));  mus_W(1:end-1,2:end-1) = mus_axy;                % Viscosities (W,E,S,N)
mus_E  = zeros(size(Vy));  mus_E(2:end  ,2:end-1) = mus_axy;
mus_S  = zeros(size(Vy));  mus_S( :     ,2:end  ) = mus;
mus_N  = zeros(size(Vy));  mus_N( :     ,1:end-1) = mus;
iVy    = NumVyG;                                                            % Eq. index for Vy (C,W,E,S,N)
iVyW   = NumVyG(1:end-1,2:end-1);
iVyE   = NumVyG(2:end  ,2:end-1);
iVyS   = NumVyG( :     ,1:end-1);
iVyN   = NumVyG( :     ,2:end  );
cVyC   = (mus_W+mus_E)/dx/dx + (mus_S+mus_N)/dy/dy;                         % Center coeff.
scVy   = max(cVyC(:));                   cVyC(:,[1,end]) = scVy;            % Scaling factor for Vy Dirichlet values
cVyW   = -mus_W(1:end-1,2:end-1)/dx/dx;                                     % West coeff.
cVyS   = -mus_S( :     ,2:end  )/dy/dy;  cVyS(:,[1,end]) = 0;               % South coeff.
Ivv    = [  iVy(:); iVyE(:); iVyN(:) ]';                                    % Triplets [I,J,V]
Jvv    = [  iVy(:); iVyW(:); iVyS(:) ]';
Vvv    = [ cVyC(:); cVyW(:); cVyS(:) ]';
%% Block VU
mus_W  = zeros(size(Vy));  mus_W( :, :     ) = musv(1:end-1, :);            % Viscosities (W,E,S,N)
mus_E  = zeros(size(Vy));  mus_E( :, :     ) = musv(2:end  , :);
mus_S  = zeros(size(Vy));  mus_S( :,2:end  ) = mus;
mus_N  = zeros(size(Vy));  mus_N( :,1:end-1) = mus;
iVy    = NumVyG( :    ,2:end-1);                                            % Eq. index for VyC
iVxSW  = NumVx(1:end-1,1:end-1);                                            % Eq. index for Vx (SW,SE,NW,NE)
iVxSE  = NumVx(2:end  ,1:end-1);
iVxNW  = NumVx(1:end-1,2:end  );
iVxNE  = NumVx(2:end  ,2:end  );
cVxSW  = (-mus_W(:,2:end-1) + mus_S(:,2:end-1))/(dx*dy);  cVxSW(1  ,:) = 0; % Coeff. for Vx (SW,SE,NW,NE)
cVxSE  = ( mus_E(:,2:end-1) - mus_S(:,2:end-1))/(dx*dy);  cVxSE(end,:) = 0;
cVxNW  = ( mus_W(:,2:end-1) - mus_N(:,2:end-1))/(dx*dy);  cVxNW(1  ,:) = 0;
cVxNE  = (-mus_E(:,2:end-1) + mus_N(:,2:end-1))/(dx*dy);  cVxNE(end,:) = 0;
Ivu    = [   iVy(:);   iVy(:);   iVy(:);   iVy(:) ]';                       % Triplets [I,J,V]
Jvu    = [ iVxSW(:); iVxSE(:); iVxNW(:); iVxNE(:) ]';
Vvu    = [ cVxSW(:); cVxSE(:); cVxNW(:); cVxNE(:) ];
cpu(4)=toc; display(['Time Build Blocks = ', num2str(cpu(4))]);
%% Assemble Blocs
tic
if SuiteSparse==1, K    = sparse2( [Iuu(:); Ivv(:); Ivu(:)], [Juu(:); Jvv(:); Jvu(:)], [Vuu(:); Vvv(:); Vvu(:)], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );
                   grad = sparse2( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny );                                      else
                   K    =  sparse( [Iuu(:); Ivv(:); Ivu(:)], [Juu(:); Jvv(:); Jvu(:)], [Vuu(:); Vvv(:); Vvu(:)], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );
                   grad =  sparse( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny );                                      end
div  = -grad';
cpu(5)=toc; display(['Time Assemble     = ', num2str(cpu(5))]);
%% BC's on K and DivV
tic
BcK                  = zeros(size(K,1),1);
BcK(NumVx(2    ,:))  = BcK(NumVx(2    ,:))  + mus(1  ,:  )'/dx/dx.*Vx_W;
BcK(NumVx(end-1,:))  = BcK(NumVx(end-1,:))  + mus(end,:  )'/dx/dx.*Vx_E;
BcK(NumVyG(:, 2   )) = BcK(NumVyG(:, 2   )) + mus(:  ,1  ) /dy/dy.*Vy_S;
BcK(NumVyG(:,end-1)) = BcK(NumVyG(:,end-1)) + mus(:  ,end) /dy/dy.*Vy_N;
BcK(ibcVxN) = BcK(ibcVxN) +  (  mus(1:end-1,end) - musv(2:end-1,end) )/dx/dy.*Vy_N(1:end-1) ;   % VyNW
BcK(ibcVxN) = BcK(ibcVxN) +  ( -mus(2:end  ,end) + musv(2:end-1,end) )/dx/dy.*Vy_N(2:end  ) ;   % VyNE
BcK(ibcVxS) = BcK(ibcVxS) +  ( -mus(1:end-1,1  ) + musv(2:end-1, 1 ) )/dx/dy.*Vy_S(1:end-1) ;   % VySW
BcK(ibcVxS) = BcK(ibcVxS) +  (  mus(2:end  ,1  ) - musv(2:end-1, 1 ) )/dx/dy.*Vy_S(2:end  ) ;   % VySE
BcK(ibcVyE) = BcK(ibcVyE) + ((  mus(end,1:end-1) - musv(end,2:end-1) )/dx/dy.*Vx_E(1:end-1)')'; % VxSE
BcK(ibcVyE) = BcK(ibcVyE) + (( -mus(end,2:end  ) + musv(end,2:end-1) )/dx/dy.*Vx_E(2:end  )')'; % VxNE
BcK(ibcVyW) = BcK(ibcVyW) + (( -mus(1  ,1:end-1) + musv(1  ,2:end-1) )/dx/dy.*Vx_W(1:end-1)')'; % VxSW
BcK(ibcVyW) = BcK(ibcVyW) + ((  mus(1  ,2:end  ) - musv(1  ,2:end-1) )/dx/dy.*Vx_W(2:end  )')'; % VxNW 
% Conforming Dirichlets
BcK([ibcVxW; ibcVxE; ibcVyS; ibcVyN]) = [Vx_W*scVx; Vx_E*scVx; Vy_S*scVy; Vy_N*scVy];
% Non-conforming Dirichlets
cNC        = [mus_axy(:,1)/dy/dy; mus_axy(:,end)/dy/dy; mus_axy(1,:)'/dx/dx; mus_axy(end,:)'/dx/dx];
d0         = spdiags(K,0);
d0(ibcNC)  = d0(ibcNC)  + 2*cNC;
BcK(ibcNC) = BcK(ibcNC) + 2*cNC.*vBcNC;
K          = spdiags(d0,0,K);
% BC on div
BcD                 = zeros(size(div,1),1);
BcD(NumPt(1  ,:  )) = BcD(NumPt(1  ,:  )) + 1/dx*Vx_W;
BcD(NumPt(end,:  )) = BcD(NumPt(end,:  )) - 1/dx*Vx_E;
BcD(NumPt( : ,1  )) = BcD(NumPt(:  ,1  )) + 1/dy*Vy_S;
BcD(NumPt( : ,end)) = BcD(NumPt(:  ,end)) - 1/dy*Vy_N;
cpu(6)=toc; display(['Time BC           = ', num2str(cpu(6))]);
%% Prepare solve
tic
K = K + K' - diag(diag(K));                                                 % Build full from tri-lower (for checking)
cpu(7)=toc;
tic
if solver==0                                                                % Build M_Stokes
    Ms  = [ K   , -div' ; ...
            div ,  PP   ];
    cpu(8)=toc; display(['Time Build Ms     = ', num2str(cpu(8))]);
else
    PPI      = spdiags(1./diag(PP),0,PP);
    Kt       = K - grad*(PPI*div);                                          % PPI*DivV = PP\DivV
    [Kc,e,s] = chol(Kt,'lower','vector');                                   % Cholesky factorization
    cpu(8)=toc; display(['Time CHOLESKY     = ', num2str(cpu(8))]);
end
tic
u   = [Vx(:) ; Vy(:)];
p   =  Pt(:);
up  = [ u ; p ];
cpu(9)=toc;
%% Solve
tic
for it=1:niter
    if solver==0,
        up  = Ms\[ BcK  ; BcD  + PP *up(NumPtG(:))];                        % Backslash
        du  = BcK -   K*up([NumVx(:); NumVyG(:)]) - grad*up(NumPtG(:));
        dp  = BcD - div*up([NumVx(:); NumVyG(:)]);
    else
        Rhs =      BcK  - grad*(PPI*BcD + p);                               % Powell-Hestenes
        if SuiteSparse==1, u(s) = cs_ltsolve(Kc,cs_lsolve(Kc,Rhs(s))); else % Powell-Hestenes
                           u(s) = Kc'\(Kc\Rhs(s));                     end  % Powell-Hestenes Matlab
        p   = p + PPI*(BcD - div*u);
        du  = BcK -   K*u - grad*p;
        dp  = BcD - div*u;
    end
    if noisy==1, fprintf('  --- iteration %d --- \n',it);
                 fprintf('   Res. |du| = %2.2e \n',norm(du)/length(du));
                 fprintf('   Res. |dp| = %2.2e \n',norm(dp)/length(dp)); end
    if norm(dp)/length(dp) < tol_lin, break, end
end%it
cpu(10)=toc; display(['Time Backsubs     = ', num2str(cpu(10))]);
tic
XPH = [u ; p];
XBS = up; 
%% Post-processing
if solver==0, up = XBS; else up = XPH; end
Pt  = reshape(up(NumPtG(:)),[nx  ,ny  ]); Pt = Pt - mean(Pt(:));
Vx  = reshape(up(NumVx(:)) ,[nx+1,ny  ]);
Vy  = reshape(up(NumVyG(:)),[nx  ,ny+1]);
u   = up([NumVx(:);NumVyG(:)]);
% Check residuals % errors
fuv = BcK -   K*u - grad*Pt(:);
fpt = BcD - div*u;
error_mom = norm(fuv)/length(fuv); display([' Error momentum        = ', num2str(error_mom)]);
error_div = norm(fpt)/length(fpt); display([' Error divergence      = ', num2str(error_div)]);
Pre = abs(Pt-Pa);
Vxe = abs(Vx-Vxa);
Vye = abs(Vy-Vya);
L1PrErr2 = sum(Pre(:)*dx*dy)./sum(ones(size(Pre(:)))*dx*dy); display([' Check Pressure vs ANALYTICAL = ', num2str(L1PrErr2)]);
L1VxErr2 = sum(Vxe(:)*dx*dy)./sum(ones(size(Vxe(:)))*dx*dy); display([' Check Vx vs ANALYTICAL       = ', num2str(L1VxErr2)]);
L1VyErr2 = sum(Vye(:)*dx*dy)./sum(ones(size(Vye(:)))*dx*dy); display([' Check Vy vs ANALYTICAL       = ', num2str(L1VyErr2)]);
figure(1),clf,colormap('jet'),set(gcf,'Color','white')
subplot(331),imagesc(flipud(Pa' )),colorbar,axis image,ylabel('Pressure'),title('Analytic')
subplot(334),imagesc(flipud(Vxa')),colorbar,axis image,ylabel('X velocity')
subplot(337),imagesc(flipud(Vya')),colorbar,axis image,ylabel('Y velocity')
subplot(332),imagesc(flipud(Pt')),axis image,colorbar,caxis([min(Pa(:)) max(Pa(:))]),title('Numeric')
subplot(335),imagesc(flipud(Vx')),axis image,colorbar
subplot(338),imagesc(flipud(Vy')),axis image,colorbar
subplot(333),imagesc(flipud(Pre')),axis image,colorbar,title('error')
subplot(336),imagesc(flipud(Vxe')),axis image,colorbar
subplot(339),imagesc(flipud(Vye')),axis image,colorbar,drawnow
% DivV and Mom X & Y Check
DivV_expl   = diff(Vx,1,1)/dx+diff(Vy,1,2)/dy;
DivV_impl   = reshape(div*u - BcD,nx,ny);
check_DivV  = norm(DivV_expl-DivV_impl); display([' Check Divergence vs explicit = ', num2str(check_DivV)]);
Vx_e        = zeros(nx+1,ny+2); Vx_e(:,2:end-1) = Vx;
Vy_e        = zeros(nx+2,ny+1); Vy_e(2:end-1,:) = Vy;
Vx_e(:,  1) = 2*Vx_S  - Vx(:,  1);
Vx_e(:,end) = 2*Vx_N  - Vx(:,end);
Vy_e(  1,:) = 2*Vy_W' - Vy(  1,:);
Vy_e(end,:) = 2*Vy_E' - Vy(end,:);
tau_xx      = 2*mus.*(diff(Vx,1,1)/dx - 1/2*DivV_expl);
tau_yy      = 2*mus.*(diff(Vy,1,2)/dy - 1/2*DivV_expl);
tau_xy      =  musv.*(diff(Vx_e,1,2)/dy + diff(Vy_e,1,1)/dx );
Res_x       = diff(-Pt + tau_xx,1,1)/dx + diff(tau_xy(2:end-1,:),1,2)/dy;
Res_y       = diff(-Pt + tau_yy,1,2)/dy + diff(tau_xy(:,2:end-1),1,1)/dx;
ResM_expl   = [Res_x(:) ; Res_y(:)];
ResM_impl   = K*u + grad*Pt(:) - BcK;
check_ResM  = abs(norm(ResM_expl)-norm(ResM_impl)); display([' Check Momentum vs explicit   = ', num2str(check_ResM)]);
figure(2),clf,colormap('jet'),set(gcf,'Color','white')
subplot(321),imagesc(flipud(DivV_expl')),colorbar,axis image,ylabel('Div'),title('Explicit Residuals')
subplot(323),imagesc(flipud(Res_x')),colorbar,axis image,ylabel('Res X')
subplot(325),imagesc(flipud(Res_y')),colorbar,axis image,ylabel('Res Y')
subplot(322),imagesc(flipud(DivV_impl')),colorbar,axis image,title('Implicit Residuals')
subplot(324),imagesc(flipud(reshape(ResM_impl(1:max(NumVx(:))),[nx+1,ny])')),colorbar,axis image
subplot(326),imagesc(flipud(reshape(ResM_impl(max(NumVx(:))+1:end),[nx,ny+1])')),colorbar,axis image,drawnow
cpu(11)=toc; display(['WALL TIME = ', num2str(sum(cpu(1:10))),' seconds (without postprocessing)']);
end
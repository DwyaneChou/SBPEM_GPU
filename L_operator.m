function [LU,LV,LZ] = L_operator(U,V,Z,dx,dy,a,Omega,g,lat_u,lat_v,lat_z,...
                                 nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                 coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)
% Prepare stagger field
h                 = sqrt(Z);
hip1(1:nx_z-1,:)  = h(2:nx_z,:);
hip1(nx_z    ,:)  = h(1,:);
him1(2:nx_z  ,:)  = h(1:nx_z-1,:);
him1(1       ,:)  = h(nx_z,:);
hjp1(:,1:ny_z-1)  = h(:,2:ny_z);
hjp1(:,ny_z    )  = mean(h(:,ny_z));
hjm1(:,2:ny_z  )  = h(:,1:ny_z-1);
hjm1(:,1       )  = mean(h(:,1));

hOnU              = 0.5*(h+him1); % h on u grid
hOnV              = 0.5*(h+hjm1); % h on v grid
hOnV(:,ny_v)      = mean(h(:,ny_z));

u                 = U./hOnU;
v                 = V./hOnV;

Zim1(2:nx_z  ,:)  = Z(1:nx_z-1,:);
Zim1(1       ,:)  = Z(nx_z,:);
Zjm1(:,2:ny_z  )  = Z(:,1:ny_z-1);
Zjm1(:,1       )  = mean(Z(:,1));

ZOnV              = Z;
ZOnV(:,ny_v)      = mean(Z(:,ny_z));
ZOnV_jm1          = Zjm1;
ZOnV_jm1(:,ny_v)  = Z(:,ny_z);

uip1(1:nx_u-1,:)    = u(2:nx_u,:);
uip1(nx_u    ,:)    = u(1,:);
ujm1(:,2:ny_u  )    = u(:,1:ny_u-1);
ujm1(:,1       )    = 0;% For polar
uip1jm1(1:nx_u-1,:) = ujm1(2:nx_u,:);
uip1jm1(nx_u    ,:) = ujm1(1,:);

uOnV                = u;
uOnV(:,ny_v)        = 0;% For polar
uOnV_ip1            = uip1;
uOnV_ip1(:,ny_v)    = 0;% For polar
uOnV_ip1jm1         = uip1jm1;
uOnV_ip1jm1(:,ny_v) = 0;% For polar
uOnV_jm1            = ujm1;
uOnV_jm1(:,ny_v)    = 0;% For polar

Uip1(1:nx_u-1,:)    = U(2:nx_u,:);
Uip1(nx_u    ,:)    = U(1,:);
Uim1(2:nx_u  ,:)    = U(1:nx_u-1,:);
Uim1(1       ,:)    = U(nx_u,:);
Ujp1(:,1:ny_u-1)    = U(:,2:ny_u);
Ujp1(:,ny_u    )    = 0;% For polar
Ujm1(:,2:ny_u  )    = U(:,1:ny_u-1);
Ujm1(:,1       )    = 0;% For polar
Uip1jm1(1:nx_u-1,:) = Ujm1(2:nx_u,:);
Uip1jm1(nx_u    ,:) = Ujm1(1,:);

Vip1(1:nx_v-1,:)  = V(2:nx_v,:);
Vip1(nx_v    ,:)  = V(1,:);
Vim1(2:nx_v  ,:)  = V(1:nx_v-1,:);
Vim1(1       ,:)  = V(nx_v    ,:);
Vjp1(:,1:ny_v-1)  = V(:,2:ny_v);
Vjp1(:,ny_v    )  = 0;% For polar
Vjm1(:,2:ny_v  )  = V(:,1:ny_v-1);
Vjm1(:,1       )  = 0;% For polar
Vim1jp1(2:nx_v,:) = Vjp1(1:nx_v-1,:);
Vim1jp1(1     ,:) = Vjp1(nx_v    ,:);

uU                  = u.*U;
uU_ip1(1:nx_u-1,:)  = uU(2:nx_u,:);
uU_ip1(nx_u    ,:)  = uU(1,:);
uU_im1(2:nx_u  ,:)  = uU(1:nx_u-1,:);
uU_im1(1       ,:)  = uU(nx_u,:);
uU_jm1(:,2:ny_u  )  = uU(:,1:ny_u-1);
uU_jm1(:,1       )  = 0;% For polar
uU_ip1jm1(:,2:ny_u) = uU_ip1(:,1:ny_u-1);
uU_ip1jm1(:,1     ) = 0;% For polar

vcos                  = v.*cos(lat_v);
vcos_im1(2:nx_v,:  )  = vcos(1:nx_v-1,:);
vcos_im1(1     ,:  )  = vcos(nx_v,:);
vcos_jp1(:,1:ny_v-1)  = vcos(:,2:ny_v);
vcos_jp1(:,ny_v    )  = 0;% For polar
vcos_im1jp1(2:nx_v,:) = vcos_jp1(1:nx_v-1,:);
vcos_im1jp1(1     ,:) = vcos_jp1(nx_v,:);

vcosOnU               = vcos(:,1:ny_u);
vcosOnU_im1           = vcos_im1(:,1:ny_u);
vcosOnU_jp1           = vcos_jp1(:,1:ny_u);
vcosOnU_im1jp1        = vcos_im1jp1(:,1:ny_u);

vVcos                 = v.*V.*cos(lat_v);
vVcos_jp1(:,1:ny_v-1) = vVcos(:,2:ny_v);
vVcos_jp1(:,ny_v    ) = 0;% For polar
vVcos_jm1(:,2:ny_v  ) = vVcos(:,1:ny_v-1);
vVcos_jm1(:,1       ) = 0;% For polar

VcosOnZ_temp          = V.*cos(lat_v);
VcosOnZ               = VcosOnZ_temp(:,1:ny_z);
VcosOnZ_jp1           = VcosOnZ_temp(:,2:ny_v);

% lat/lon coef.
lat_u_jm1(:,2:ny_u)   = lat_u(:,1:ny_u-1);
lat_u_jm1(:,1     )   = 0;% For polar

sinLatU               = sin(lat_u);
sinLatU_jm1           = sin(lat_u_jm1);
tanLatU               = tan(lat_u);
tanLatU_jm1           = tan(lat_u_jm1);

lat_v_jp1(:,1:ny_v-1) = lat_v(:,2:ny_v);
lat_v_jp1(:,ny_v    ) = lat_v(:,ny_v);

cosLatV               = cos(lat_v);
cosLatVOnU            = cosLatV(:,1:ny_u);
cosLatV_jp1           = cos(lat_v_jp1);
cosLatVOnU_jp1        = cosLatV_jp1(:,1:ny_u);

cosLatZ               = cos(lat_z);

% Advection
ADV_U_x  = u      .*(Uip1-Uim1) + uU_ip1    - uU_im1;

ADV_U_y  = (vcosOnU_jp1 + vcosOnU_im1jp1).*Ujp1 - (vcosOnU + vcosOnU_im1).*Ujm1;

ADV_V_x  = (uOnV_ip1    + uOnV_ip1jm1   ).*Vip1 - (uOnV    + uOnV_jm1   ).*Vim1;

ADV_V_y  = vcos   .*(Vjp1-Vjm1) + vVcos_jp1 - vVcos_jm1;

% Flux
FLUX_Z_x = (hip1+h).*Uip1        - (h+him1).*U;
FLUX_Z_y = (hjp1+h).*VcosOnZ_jp1 - (h+hjm1).*VcosOnZ;

% Pressure Gradient Force
PGF_U_temp        = hOnU.*(Z-Zim1);
PGF_U             = 4.0*coefU_x.*PGF_U_temp;

PGF_V_temp        = hOnV.*(ZOnV-ZOnV_jm1);
PGF_V             = PGF_V_temp/(a*dy);
PGF_V(:,ny_v)     = 0;
PGF_V(:,1   )     = 0;

% Colioris Force
C1                 = cosLatVOnU_jp1./cosLatZ;
C2                 = cosLatVOnU    ./cosLatZ;
VOnU_part1_temp    = Vjp1+Vim1jp1;
VOnU_part1         = VOnU_part1_temp(:,1:ny_u).*C1;
VOnU_part2_temp    = V+Vim1;
VOnU_part2         = VOnU_part2_temp(:,1:ny_u).*C2;
VOnU               = 0.25*(VOnU_part1+VOnU_part2);

fv                 = 2.0*Omega*sinLatU.*VOnU;

fu                 = 0.5*Omega*((Uip1+U).*sinLatU+(Uip1jm1+Ujm1).*sinLatU_jm1); % 2*Omega*sin(theta)U/4
fu(:,ny_v)         = 0;
fu(:,1   )         = 0;

% Curvature Term
CV                 = u/a.*tanLatU.*VOnU; % u/a*tan(theta)*V_On_U

CU_part1           = (uU_ip1   +uU    ).*tanLatU;
CU_part2           = (uU_ip1jm1+uU_jm1).*tanLatU_jm1;
CU                 = 0.25*(CU_part1+CU_part2)/a;
CU(:,1   )         = 0;
CU(:,ny_v)         = 0;

% Construct L operator
LU1         = ADV_U_x;
LU2         = ADV_U_y;
LU          = coefU_x.*LU1+coefU_y.*LU2+PGF_U-fv-CV;

LV1         = ADV_V_x;
LV2         = ADV_V_y;
LV1(:,1   ) = 0;
LV1(:,ny_v) = 0;
LV2(:,1   ) = 0;
LV2(:,ny_v) = 0;
LV          = coefV_x.*LV1+coefV_y.*LV2+PGF_V+fu+CU;

LZ1         = coefZ_x.*FLUX_Z_x;
LZ2         = coefZ_y.*FLUX_Z_y;
LZ          = LZ1+LZ2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Symmetry Check %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% LUU         = LU.*U;
% LVV         = LV.*V;
% LZZ         = LZ.*Z;
% 
% figure
% pcolor(LUU)
% shading interp
% 
% figure
% pcolor(LVV)
% shading interp
% 
% figure
% pcolor(LZZ)
% shading interp

% CFUU = fv.*U;
% CFVV = fu.*V;
% CTUU = CV.*U;
% CTVV = CU.*V;
% 
% LU1_U = sum(sum(LU1.*U));
% LU2_U = sum(sum(LU2.*U));
% LV1_V = sum(sum(LV1.*V));
% LV2_V = sum(sum(LV2.*V));
% PGUU  = sum(sum(PGF_U_temp.*U));
% LZZx  = sum(sum(FLUX_Z_x.*Z));
% PGVV  = sum(sum(PGF_V_temp.*V.*cos(lat_v)));
% LZZy  = sum(sum(FLUX_Z_y.*Z));
% fVU   = sum(sum(CFUU.*cos(lat_u)));
% fUV   = sum(sum(CFVV.*cos(lat_v)));
% CUV   = sum(sum(CTUU.*cos(lat_u)));
% CVU   = sum(sum(CTVV.*cos(lat_v)));
% 
% disp(['(LU1,U)    = ',num2str(LU1_U)])
% disp(['(LU2,U)    = ',num2str(LU2_U)])
% disp(['(LV1,V)    = ',num2str(LV1_V)])
% disp(['(LV2,V)    = ',num2str(LV2_V)])
% disp(['(PGUU+LZZ) = ',num2str(PGUU+0.5*LZZx)]);
% disp(['(PGVV+LZZ) = ',num2str(PGVV+0.5*LZZy)]);
% 
% disp(['(fVU+fUV) = ',num2str(-fVU+fUV)]);
% disp(['(CUV+CVU) = ',num2str(-CUV+CVU)]);
% disp('                                  ');
% ;
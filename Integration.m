function [tau_n,U_np1,V_np1,Z_np1] = Integration(time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                                       lat_u,lat_v,lat_z,...
                                                       nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                                       coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)

[LU,LV,LZ,BU,BV,BZ] = B_operator(time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                 lat_u,lat_v,lat_z,...
                                 nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                 coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
        
% U
cosU     = cos(lat_u);
LU_norm2 = sum(sum(LU.*LU.*cosU));
BUU      = sum(sum(BU.*U .*cosU));
BULU     = sum(sum(BU.*LU.*cosU));
BU_norm2 = sum(sum(BU.*BU.*cosU));

% V
cosV     = cos(lat_v);
LV_norm2 = sum(sum(LV.*LV.*cosV));
BVV      = sum(sum(BV.*V .*cosV));
BVLV     = sum(sum(BV.*LV.*cosV));
BV_norm2 = sum(sum(BV.*BV.*cosV));

% Z
cosZ     = cos(lat_z);
LZ_norm2 = sum(sum(LZ.*LZ.*cosZ));
BZZ      = sum(sum(BZ.*Z .*cosZ));
BZLZ     = sum(sum(BZ.*LZ.*cosZ));
BZ_norm2 = sum(sum(BZ.*BZ.*cosZ));

% F
LF_norm2 = LU_norm2 + LV_norm2 + LZ_norm2;
BFF      = BUU      + BVV      + BZZ;
BFLF     = BULU     + BVLV     + BZLZ;
BF_norm2 = BU_norm2 + BV_norm2 + BZ_norm2;

K1       = -LF_norm2/BFF;
K2       = BFLF/BFF;
K3       = -BF_norm2/BFF;

eps      = gpuArray(0.5);
tau1     = 2.0 - K1 / eps;
tau2     = K2 + sqrt( K2^2 + eps * tau1 * K3 );
tau_n    = tau1/tau2;

U_np1    = U - LU * tau_n + eps * tau_n.^2 * BU;
V_np1    = V - LV * tau_n + eps * tau_n.^2 * BV;
Z_np1    = Z - LZ * tau_n + eps * tau_n.^2 * BZ;
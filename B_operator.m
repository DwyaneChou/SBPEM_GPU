function [LU1,LV1,LZ1,BU,BV,BZ] = B_operator(time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                             lat_u,lat_v,lat_z,...
                                             nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                             coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)

U0 = U;
V0 = V;
Z0 = Z;

% K1
[LU,LV,LZ] = L_operator(U,V,Z,dlambda,dtheta,a,Omega,g,...
                        lat_u,lat_v,lat_z,...
                        nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                        coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
LU1  = LU;
LV1  = LV;
LZ1  = LZ;

U_k1 = -LU;
V_k1 = -LV;
Z_k1 = -LZ;

% % Inverse symmetry check
% LUU  = sum(sum(LU.*U.*cos(lat_u)));
% LVV  = sum(sum(LV.*V.*cos(lat_v)));
% LZZ  = sum(sum(LZ.*Z.*cos(lat_z)));
% disp(num2str(LUU+LVV+LZZ))

% K2
U    = U0 + 0.5*time_step*U_k1;
V    = V0 + 0.5*time_step*V_k1;
Z    = Z0 + 0.5*time_step*Z_k1;

[LU,LV,LZ] = L_operator(U,V,Z,dlambda,dtheta,a,Omega,g,...
                        lat_u,lat_v,lat_z,...
                        nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                        coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
U_k2 = -LU;
V_k2 = -LV;
Z_k2 = -LZ;

% % Inverse symmetry check
% LUU  = sum(sum(LU.*U.*cos(lat_u)));
% LVV  = sum(sum(LV.*V.*cos(lat_v)));
% LZZ  = sum(sum(LZ.*Z.*cos(lat_z)));
% disp(num2str(LUU+LVV+LZZ))

% K3
U    = U0 + 0.5*time_step*U_k2;
V    = V0 + 0.5*time_step*V_k2;
Z    = Z0 + 0.5*time_step*Z_k2;

[LU,LV,LZ] = L_operator(U,V,Z,dlambda,dtheta,a,Omega,g,...
                        lat_u,lat_v,lat_z,...
                        nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                        coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
U_k3 = -LU;
V_k3 = -LV;
Z_k3 = -LZ;

% % Inverse symmetry check
% LUU  = sum(sum(LU.*U.*cos(lat_u)));
% LVV  = sum(sum(LV.*V.*cos(lat_v)));
% LZZ  = sum(sum(LZ.*Z.*cos(lat_z)));
% disp(num2str(LUU+LVV+LZZ))

% K4
U    = U0 + time_step*U_k3;
V    = V0 + time_step*V_k3;
Z    = Z0 + time_step*Z_k3;

[LU,LV,LZ] = L_operator(U,V,Z,dlambda,dtheta,a,Omega,g,...
                        lat_u,lat_v,lat_z,...
                        nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                        coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
U_k4  = -LU;
V_k4  = -LV;
Z_k4  = -LZ;

% % Inverse symmetry check
% LUU  = sum(sum(LU.*U.*cos(lat_u)));
% LVV  = sum(sum(LV.*V.*cos(lat_v)));
% LZZ  = sum(sum(LZ.*Z.*cos(lat_z)));
% disp(num2str(LUU+LVV+LZZ))

phi4_U = (U_k1 + 2.0*U_k2 + 2.0*U_k3 + U_k4)/6.0;
phi4_V = (V_k1 + 2.0*V_k2 + 2.0*V_k3 + V_k4)/6.0;
phi4_Z = (Z_k1 + 2.0*Z_k2 + 2.0*Z_k3 + Z_k4)/6.0;

BU     = 2.0 / time_step * (phi4_U - U_k1);
BV     = 2.0 / time_step * (phi4_V - V_k1);
BZ     = 2.0 / time_step * (phi4_Z - Z_k1);
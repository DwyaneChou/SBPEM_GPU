function [u,v,Z] = Haurwitz(a,Omega,g,lon_u,lat_u,lon_v,lat_v,lon_z,lat_z)
% Set Constants
omega = 7.848*10^-6;
K     = omega;
R     = 4.0;
h0    = 8.0*10^3;

% Initial fields with Rossby-Haurwitz Wave, reference:
% "A Standard Test Set for Numerical Approximations to the Shallow
% Water Equations in Spherical Geometry"
u1  = a*omega*cos(lat_u);
u21 = a*K*cos(lat_u).^(R-1.0);
u22 = R.*sin(lat_u).^2-cos(lat_u).^2;
u23 = cos(R*lon_u);
u   = u1+u21.*u22.*u23;

v = -a*K*R*cos(lat_v).^(R-1).*sin(lat_v).*sin(R*lon_v);

A1  = omega*0.5*(2*Omega+omega)*cos(lat_z).^2;
Ac  = 0.25*K^2*cos(lat_z).^(2.0*R);
A21 = (R+1.0).*cos(lat_z).^2;
A22 = 2.0*R^2-R-2;
A23 = 2.0*R^2.*cos(lat_z).^-2;
A   = A1+Ac.*(A21+A22-A23);

Bc  = 2.*(Omega+omega)*K/((R+1)*(R+2)).*cos(lat_z).^R;
B1  = R^2+2.0*R+2.0;
B2  = (R+1.0)^2.*cos(lat_z).^2;
B   = Bc.*(B1-B2);

Cc  = 0.25*K^2.*cos(lat_z).^(2.0*R);
C1  = (R+1.0)*cos(lat_z).^2;
C2  = R+2.0;
C   = Cc.*(C1-C2);

Z  = g*h0+a^2*A+a^2*B.*cos(R*lon_z)+a^2*C.*cos(2.0*R*lon_z);
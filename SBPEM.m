% Spherical Barotropic Primitive Eqation Model
% Written by Zhou Lilong Jun 22,2018

clc
clear

time_start = clock;

% Define Constants
Omega = gpuArray(7.292*10^-5);
a     = gpuArray(6371229.0);
g     = gpuArray(9.8);

% Define the grid resolution
dx = 0.5; % Degree
dy = 0.5; % Degree

% Define time(seconds)
max_time_step = cosd(90-dy)*a/180.0*pi/1000*6;
time_step     = gpuArray(1);
run_time      = gpuArray(33*24*3600);

% Define output
history_interval = gpuArray(3600);
output_precision = 'NC_FLOAT';

% Generate C-grid on sphere, z represents the geopotential height
longitude_z = 0:dx:360-dx;
latitude_z  = -90+0.5*dx:dx:90-0.5*dx;

longitude_u = 0+0.5*dx:dx:360+0.5*dx-dx;
latitude_u  = -90+0.5*dx:dx:90-0.5*dx;

longitude_v = 0:dx:360-dx;
latitude_v  = -90:dx:90;

[lat_u,lon_u] = meshgrid(latitude_u,longitude_u);
[lat_v,lon_v] = meshgrid(latitude_v,longitude_v);
[lat_z,lon_z] = meshgrid(latitude_z,longitude_z);

% Convert longitude/latitude from degree to radian
d2r   = pi/180.0;
lon_u = gpuArray(lon_u*d2r);
lat_u = gpuArray(lat_u*d2r);
lon_v = gpuArray(lon_v*d2r);
lat_v = gpuArray(lat_v*d2r);
lon_z = gpuArray(lon_z*d2r);
lat_z = gpuArray(lat_z*d2r);

dlambda = gpuArray(dx*d2r);
dtheta  = gpuArray(dy*d2r);

% % Plot V grid
% [x,y,z] = sph2cart(lon_v,lat_v,a);
% surf(x,y,z);

% Get grid size
nx_u = gpuArray(size(lon_u,1));
ny_u = gpuArray(size(lon_u,2));
nx_v = gpuArray(size(lon_v,1));
ny_v = gpuArray(size(lon_v,2));
nx_z = gpuArray(size(lon_z,1));
ny_z = gpuArray(size(lon_z,2));

% Initial fields with Rossby-Haurwitz Wave
[u,v,Z] = Haurwitz(a,Omega,g,lon_u,lat_u,lon_v,lat_v,lon_z,lat_z);
u       = gpuArray(u);
v       = gpuArray(v);
Z       = gpuArray(Z);

% IAP transformation
h                 = sqrt(Z);
him1(2:nx_z  ,:)  = h(1:nx_z-1,:);
him1(1       ,:)  = h(nx_z,:);
hjm1(:,2:ny_z  )  = h(:,1:ny_z-1);
hjm1(:,1       )  = mean(h(:,1));

hOnU              = 0.5*(h+him1); % h on u grid
hOnV              = 0.5*(h+hjm1); % h on v grid
hOnV(:,ny_v)      = mean(h(:,ny_z));

U                 = hOnU.*u;
V                 = hOnV.*v;

total_energy0     = sum(sum(U.*U.*cos(lat_u)))+sum(sum(V.*V.*cos(lat_v)))+sum(sum(Z.*Z.*cos(lat_z)));
total_mass0       = sum(sum(Z.*cos(lat_z)));

% Compute the coefficient for L operator
coefU_x = 0.25./(a*cos(lat_u)*dlambda); %coefU_x = 1/(2*a*cos(theta_u)*2dx)
coefU_y = 0.25./(a*cos(lat_u)*dtheta ); %coefU_y = 1/(2*a*cos(theta_u)*2dy)
coefV_x = 0.25./(a*cos(lat_v)*dlambda); %coefV_x = 1/(2*a*cos(theta_v)*2dx)
coefV_y = 0.25./(a*cos(lat_v)*dtheta ); %coefV_x = 1/(2*a*cos(theta_v)*2dy)
coefZ_x = 0.5 ./(a*cos(lat_z)*dlambda); %coefZ_x = 1/(  a*cos(theta_z)*2dx)
coefZ_y = 0.5 ./(a*cos(lat_z)*dtheta ); %coefZ_x = 1/(  a*cos(theta_z)*2dy)

% Output the initial status
int_step_num = ceil(run_time/time_step);
output_count = 0;
output_num   = ceil(run_time/history_interval)+1;

output_netCDF(output_num,output_count,output_precision,...
              U,V,Z,dlambda,dtheta,a,Omega,g,...
              lon_u,lon_v,lon_z,lat_u,lat_v,lat_z,...
              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z)
          
ti_start = clock;
for nt = 1:int_step_num
    [time_step_n,U,V,Z] = Integration(time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                      lat_u,lat_v,lat_z,...
                                      nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                      coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
    
    integral_time = nt*time_step;
    
    % Output
	if rem(integral_time,history_interval)==0 && integral_time>=history_interval
        ti_end = clock;
        
        output_count = output_count+1;
        output_netCDF(output_num,output_count,output_precision,...
                      U,V,Z,dlambda,dtheta,a,Omega,g,...
                      lon_u,lon_v,lon_z,lat_u,lat_v,lat_z,...
                      nx_u,ny_u,nx_v,ny_v,nx_z,ny_z)
                  
        total_energy       = sum(sum(U.*U.*cos(lat_u)))+sum(sum(V.*V.*cos(lat_v)))+sum(sum(Z.*Z.*cos(lat_z)));
        total_mass         = sum(sum(Z.*cos(lat_z)));
        energy_change_rate = (total_energy-total_energy0)/total_energy0; %ECR
        mass_change_rate   = (total_mass-total_mass0)/total_mass0;       %MCR
        
        disp(['Output hour   = ',num2str(integral_time/3600),'/',num2str(output_num-1)]);
        disp(['tau_n         = ',num2str(time_step_n)])
        disp(['Total Energy  = ',num2str(total_energy)]);
        disp(['Total Mass    = ',num2str(total_mass)]);
        disp(['ECR           = ',num2str(energy_change_rate)]);
        disp(['MCR           = ',num2str(mass_change_rate)]);
        disp(['Integral time = ',num2str(etime(ti_end,ti_start))]);
        disp('                                     ');
        
        ti_start = clock;
	end
end

time_end = clock;
disp(['It took ',num2str(etime(time_end,time_start)),' seconds to run SBPEM'])
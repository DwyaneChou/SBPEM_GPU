function output_netCDF(int_step_num,output_num,output_precision,...
                       U,V,Z,dlambda,dtheta,a,Omega,g,...
                       lon_u,lon_v,lon_z,lat_u,lat_v,lat_z,...
                       nx_u,ny_u,nx_v,ny_v,nx_z,ny_z)

f_out   = 'output.nc';
    
r2d     = 180.0/pi;
lon_u   = lon_u*r2d;
lat_u   = lat_u*r2d;
lon_v   = lon_v*r2d;
lat_v   = lat_v*r2d;
lon_z   = lon_z*r2d;
lat_z   = lat_z*r2d;
dlambda = dlambda*r2d;
dtheta  = dtheta*r2d;

% Inverse IAP transformation
h                 = sqrt(Z);
him1(2:nx_z  ,:)  = h(1:nx_z-1,:);
him1(1       ,:)  = h(nx_z,:);
hjm1(:,2:ny_z  )  = h(:,1:ny_z-1);
hjm1(:,1       )  = mean(h(:,1));

hOnU              = 0.5*(h+him1); % h on u grid
hOnV              = 0.5*(h+hjm1); % h on v grid
hOnV(:,ny_v)      = mean(h(:,ny_z));

u                 = U./hOnU;
v                 = V./hOnV;

% Write Data into netCDF
west_east        = gather(nx_u);
south_north      = gather(ny_u);
south_north_stag = gather(ny_v);

if output_num==0
    mode           = netcdf.getConstant('NETCDF4');
    mode           = bitor(mode,netcdf.getConstant('CLOBBER'));
    ncid           = netcdf.create(f_out,mode);
    disp(['ncid = ',num2str(ncid)])
    
    % Define Dimensions
    time_dimID             = netcdf.defDim(ncid,'time'            ,gather(int_step_num));
    west_east_dimID        = netcdf.defDim(ncid,'west_east'       ,gather(west_east));
    south_north_dimID      = netcdf.defDim(ncid,'south_north'     ,gather(south_north));
    south_north_stag_dimID = netcdf.defDim(ncid,'south_north_stag',gather(south_north_stag));
    
    % Define Attribute
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model_Source' ,'SBPEM written by Zhou Lilong')
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dlambda'      ,gather(dlambda))
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dtheta'       ,gather(dtheta))
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'earth_radius' ,gather(a))
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Omega'        ,gather(Omega))
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'g'            ,gather(g))
    
    % Define Variables
    XLONG_U_id = netcdf.defVar(ncid,'XLONG_U',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLONG_U_id,'units','degrees longitude');
    netcdf.putAtt(ncid,XLONG_U_id,'description','longitude for u');
    
    XLAT_U_id = netcdf.defVar(ncid,'XLAT_U',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLAT_U_id,'units','degrees latitude');
    netcdf.putAtt(ncid,XLAT_U_id,'description','latitude for u');
    
    XLONG_V_id = netcdf.defVar(ncid,'XLONG_V',output_precision,[west_east_dimID,south_north_stag_dimID]);
    netcdf.putAtt(ncid,XLONG_V_id,'units','degrees longitude');
    netcdf.putAtt(ncid,XLONG_V_id,'description','longitude for v');
    
    XLAT_V_id = netcdf.defVar(ncid,'XLAT_V',output_precision,[west_east_dimID,south_north_stag_dimID]);
    netcdf.putAtt(ncid,XLAT_V_id,'units','degrees latitude');
    netcdf.putAtt(ncid,XLAT_V_id,'description','latitude for v');
    
    XLONG_M_id = netcdf.defVar(ncid,'XLONG_M',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLONG_M_id,'units','degrees longitude');
    netcdf.putAtt(ncid,XLONG_M_id,'description','longitude for z');
    
    XLAT_M_id = netcdf.defVar(ncid,'XLAT_M',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLAT_M_id,'units','degrees latitude');
    netcdf.putAtt(ncid,XLAT_M_id,'description','latitude for z');
    
    U_id = netcdf.defVar(ncid,'U',output_precision,[west_east_dimID,south_north_dimID,time_dimID]);
    netcdf.putAtt(ncid,U_id,'units','m/s');
    netcdf.putAtt(ncid,U_id,'description','u wind component');
    
    V_id = netcdf.defVar(ncid,'V',output_precision,[west_east_dimID,south_north_stag_dimID,time_dimID]);
    netcdf.putAtt(ncid,V_id,'units','m/s');
    netcdf.putAtt(ncid,V_id,'description','v wind component');
    
    Z_id = netcdf.defVar(ncid,'Z',output_precision,[west_east_dimID,south_north_dimID,time_dimID]);
    netcdf.putAtt(ncid,Z_id,'units','m^2/s^2');
    netcdf.putAtt(ncid,Z_id,'description','geopotential height');
    
    % Put Variables
    netcdf.putVar(ncid, XLONG_U_id ,gather(lon_u));
    netcdf.putVar(ncid, XLAT_U_id  ,gather(lat_u));
    netcdf.putVar(ncid, XLONG_V_id ,gather(lon_v));
    netcdf.putVar(ncid, XLAT_V_id  ,gather(lat_v));
    netcdf.putVar(ncid, XLONG_M_id ,gather(lon_z));
    netcdf.putVar(ncid, XLAT_M_id ,gather(lat_z));
    netcdf.putVar(ncid, U_id     ,[0,0,0],[west_east,south_north,1]     ,gather(u));
    netcdf.putVar(ncid, V_id     ,[0,0,0],[west_east,south_north_stag,1],gather(v));
    netcdf.putVar(ncid, Z_id     ,[0,0,0],[west_east,south_north,1]     ,gather(Z));
    
    netcdf.close(ncid)
    
else
    ncid = netcdf.open(f_out,'WRITE');
    
    U_id = netcdf.inqVarID(ncid,'U');
    V_id = netcdf.inqVarID(ncid,'V');
    Z_id = netcdf.inqVarID(ncid,'Z');
    
    netcdf.reDef(ncid)
    
    netcdf.putVar(ncid, U_id     ,[0,0,output_num],[west_east,south_north,1]     ,gather(u));
    netcdf.putVar(ncid, V_id     ,[0,0,output_num],[west_east,south_north_stag,1],gather(v));
    netcdf.putVar(ncid, Z_id     ,[0,0,output_num],[west_east,south_north,1]       ,gather(Z));
    
    netcdf.close(ncid)
end


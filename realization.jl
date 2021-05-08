function realization(ngrid,H,Br,m,y,bx,by,pd,q)

# computes the velocity and acceleration on a lat,long grid
# at a series of snapshots.
#
# input
#
# ngrid -  number of grid points in colatitude
# H    - layer thickness (km)
# Np   - dimensionless buoyancy frequency (N/Omega)
# Br   - radial magnetic field  (mT)
# m    - angular order m
# y    - meridional coordinate
# bx,by - x and y components of magnetic perturbation
# pd    - period  (years)
# q     - quality factor
#
# output
#
# tvec      - array to times to evaluate velocity
# lat,long  -  array of latitude and longitude
# vy,vx     -  theta and phi components of velocity at CMB
# ay,ax    -  theta and phi components of acceleration


# array of times
dt = 0.5;
tvec = (0.0 : dt : abs(pd));

# latitude and longitude grid
#ngrid = 100;
dlat = pi / ngrid;
dlong = pi / ngrid;
lat = (-pi/2+dlat/2 : dlat : pi/2-dlat/2);
long = (dlong/2 : dlong : 2 * pi-dlong/2);
colat = (dlat/2 : dlat : pi-dlat/2);  # reverse order

# allocate memory for velocity and acceleration
vx = Array{Float64,3}(undef,length(tvec),ngrid,2*ngrid);
vy = Array{Float64,3}(undef,length(tvec),ngrid,2*ngrid);
ax = Array{Float64,3}(undef,length(tvec),ngrid,2*ngrid);
ay = Array{Float64,3}(undef,length(tvec),ngrid,2*ngrid);

# loop over times
for k = 1 : length(tvec)

   tstep = tvec[k];
   vxt,vyt,axt,ayt = velocity(H,Br,m,y,lat,long,bx,by,pd,q,tstep);
   vx[k,:,:]=vxt;
   vy[k,:,:]=vyt;
   ax[k,:,:]=axt;
   ay[k,:,:]=ayt;

end

# collect into arrays
t = collect(tvec);
clt = collect(colat);
lg = collect(long);

return t,clt,lg,vx,vy,ax,ay

end

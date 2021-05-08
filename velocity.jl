function velocity(H,Br,m,y,lat,long,bx,by,pd,q,t)

# computes the velocity field on a lat,long grid
# at time t using a mode from pWaves().
#
# input
#
# H    - layer thickness (km)
# Np   - dimensionless buoyancy frequency (N/Omega)
# Br   - radial magnetic field  (mT)
# m    - angular order
# y    - meridional coordinate
# lat  - grid in latitude   (radian)
# long - grid in longitude  (radian)
# bx,by - x and y components of magnetic perturbation
# pd    - period  (years)
# q     - quality factor
# t     - time (in years)
#
# output
#
# vy,vx     -  theta and phi components of velocity at CMB
# ay,ax    -  theta and phi components of acceleration

# reconstruct frequency of mode (in seconds)
sectoyear = 365.25 * 24 * 60 * 60;
tsec = t * sectoyear;
omega0 = 2 * pi / (pd*sectoyear);
omega = omega0 - abs(0.5*omega0/q)im;

# set model parameters (compare with ModelParameters.jl)
B = Br * 1.0e-3;
mu = 4 * pi * 1.0e-7;      # permeability
sigma = 1.0e6;              # conductivity (S/m)
eta = 1.0 / (mu * sigma);   # diffusivity (m^2/s)
kz = pi / (H * 1000.0);      # radial wavenumber
chi = 1.0 + (eta * kz^2/omega)im;


# add southern hemisphere
nlat = length(lat);
nsol = length(y);
yw = [-reverse(y);y];
latw = asin.(yw);
bxw = [zeros(nsol);bx];
byw = [zeros(nsol);by];

#  interpolate wave on to colat,long grid
nodes = (latw,);
itp = Interpolations.interpolate(nodes,bxw,Gridded(Linear()));
bxgrid = itp(lat);
itp= Interpolations.interpolate(nodes,byw,Gridded(Linear()));
bygrid = itp(lat);

# evaluate velocity and acceleration on grid
vx = Array{Float64,2}(undef,length(lat),length(long));
vy = Array{Float64,2}(undef,length(lat),length(long));
ax = Array{Float64,2}(undef,length(lat),length(long));
ay = Array{Float64,2}(undef,length(lat),length(long));

# store in reverse order in latitude
for j = 1 : length(lat)
    for k = 1 : length(long)

       phase = exp((m*long[k] - omega*tsec)im);
       vx[nlat-j+1,k] = real(((omega*chi*bxgrid[j]*phase)/(kz*B))im);
       vy[nlat-j+1,k] = real(((omega*chi*bygrid[j]*phase)/(kz*B))im);
       ax[nlat-j+1,k] = real((omega^2*chi*bxgrid[j]*phase)/(kz*B));
       ay[nlat-j+1,k] = real((omega^2*chi*bygrid[j]*phase)/(kz*B));

    end
end

# normalize and convert to km/yr and km/yr^2
sectoyear = 365.25 * 24 * 60 * 60;
vymax = maximum(vy) * sectoyear;
aymax = maximum(ay) * sectoyear^2;

vx *= sectoyear/vymax;
vy *= sectoyear/vymax;
ax *= sectoyear^2/vymax;
ay *= sectoyear^2/vymax;


return vx,vy,ax,ay

end

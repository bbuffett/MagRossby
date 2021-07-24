function pWaves(ngrid,H,Np,Br,m,T)

# Evaluates the frequency and spatial structure of
# polar magnetic Rossby waves in a spherical geometry.
# The wave equation is formulated in terms of the
# meridional component of the magnetic field.
# The meridional coordinate is y = cos(colat), although
# a second grid array is stored as latitude for plotting.
# The radial magnetic field is assumed to be constant
# over the surface of the CMB. The same governing equation
# describe equatorial MAC waves in a spherical geometry. The
# two types of waves are distinguished by the sign of the frequency.
#
# input
#
# ngrid - grid used to represent the waves
# H    - layer thickness (km)
# Np   - dimensionless buoyancy frequency (N/Omega)
# Br   - radial magnetic field  (mT)
# m    - angular order m
# T    - initial estimate of period  (years, -'ve' for westward)
#
# output
#
# lat  -  array of latitude for plotting eigenfunction by(y)
# y    -  array of coordinates "y"
# byp   -  eigenfunction by(y)  (modified according to Eq. 13)
# omega - frequency of wave
# pd    - period of wave (years)
# q     - quality factor

# initial estimate of frequency
sectoyear = 365.25 * 24 * 60 * 60;
omega0 = 2 * pi / (T*sectoyear) + 0.0im;
B = Br * 1.0e-3;    # convert to Tesla

# set up grid
y,lat = mgrid(ngrid);

# solve eigenvalue problem for byp
omega,byp = modes_y(y,H,Np,B,m,omega0);

#  compute b_x
bxp = modes_x(y,byp,H,Np,B,m,omega);

# transform back to original variables
ynew,bx,by = transform(y,bxp,byp);

# evaluate period and quality factor of wave
pd = 2 * pi / (sectoyear * real(omega));
q = abs(0.5*real(omega)/imag(omega));

return pd,q,bx,by,ynew

end

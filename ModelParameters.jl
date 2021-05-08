function  ModelParameters(H,Np,B,m,omega)

# This function evaluates the constants needed in
# the MAC-wave model

# input
#   H      - layer thickness (km)
#   Np     - dimensionless buoyancy frequency (N/Omega)
#   B      - radial magnetic field (T)
#   m      - angular order
#   omega  - nominal wave frequency
#   omega0 - reference frequency for shift - invert

# output
#   C  -  Coriolis parameter in wave equation
#   M  -  magnetic parameter in wave equation

#
# some preliminary constants
#Br = 0.6e-3;               # radial magnetic field (T)
mu = 4 * pi * 1.0e-7;      # permeability
sigma = 1.0e6;              # conductivity (S/m)
eta = 1.0 / (mu * sigma);   # diffusivity (m^2/s)
Omega = 0.7292e-4;          # rotation rate
rho = 1.0e4;                # density  (kg/m^3)
N = Np * Omega;             # stratification
R = 3.48e6;                 # radius of core

# intermediate results
Va = B / sqrt(mu * rho);
kz =  pi / (1000 * H);

# evaluate diffuion factor using estimate of omega
 chi = 1  +  eta * kz^2 * im / (omega);

 # model constants
 C = 2*Omega*omega*(kz*R/N)^2;
 M = Va^2 * kz^2 * (kz*R/N)^2 / chi;

return C, M


end

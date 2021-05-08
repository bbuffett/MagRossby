function modes_x(y,byp,H,Np,B,m,omega)

# Computes the longitudinal component of the magnetic
# perturbation using a solution for the meridional
# component
#
# input
#  y     - grid points in meridional coordinate
#  byp   - meridional magnetic perturbation
#  H     - thickness of layer (km)
#  Np    - dimensionless buoyancy frequency  (N/Omega)
#  B     - radial magnetic field  (T)
#  m     - angular order
#  omega0 - initial estimate of frequency


# a few constants   (must be consistent with ModelParameters)
Omega = 0.7292e-4;
R = 3.48e6;
N = Np * Omega;
k = pi/(1000.0*H);

# compute model parameters
C,M = ModelParameters(H,Np,B,m,omega);

# evaluate derivative matrices
d1,d2 = derivatives(y);

#   right-hand side of Eq. (17) for bxp
n = length(y);
y2 = y.*y;
yc2 = -(y2 .- 1.0);
a1 = -(m^2*ones(n) + M * yc2);

#   left-hand side
a2 = -((C)im * y .* byp);
a3 = (m * d1*byp)im;

bxp = (a2 + a3)./a1;


return bxp

end

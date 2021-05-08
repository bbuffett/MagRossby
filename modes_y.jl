function modes_y(y,H,Np,B,m,omega0)

# This function constructs and solves the eigenvalue problem
# using NonlinearEigenproblems package
#
# input
#  y     - grid points in meridional coordinate
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
C0,M = ModelParameters(H,Np,B,m,omega0);

# evaluate coefficients of wave equation
a0,a1,a2 = wave_equation(m,M,y);

# eigenvalue problem
nep = PEP([a0,a1,a2]);
d,v = polyeig(nep);

# which eigenvalue?
index = findmin(abs.(real(d) .- C0))[2];
println("mode number = ",index)

# unpack frequency
Cnew = d[index];
by = v[:,index];
wnew = Cnew * N^2 / (2*Omega*(k*R)^2);

return wnew,by

end

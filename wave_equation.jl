function wave_equation(m,M,y)

# evaluates the discrete matrix equation for the waves
#
# a1 b'' + a2 b' + (a3 + a4 lambda + a5 I lambda^2)b = 0
#
# where a1-a2 are tridiagonal matrices resulting from derivatives
# and a3-a5 are diagonal matrices. The eigenvalue is lambda = C
# (see eq. 37) B&M 2019. We recover the value of frequency from
# the definition of C
#
# input
#  m -  angular order of wave
#  C -  Coriolis parameter  (dimensionless)
#  M -  magnetic parameter  (dimensionless)
#  y -  meridional coordinate
#
# output
# A0,A1,A2 - matrices used in nonlinear eigenvalue problem
#
# Form of eigenvalue problem
# (A0 + A1 lambda + A2 lambda^2 ) b = 0
#

# size of matrix system
n = length(y);

# evaluate first and second derivatives
d1,d2 = derivatives(y);

# assemble matrix a1 (second derivative)
y2 = y.*y;
y2c = -(y2.-1);
a1 = Diagonal(y2c) * d2;

# assemble matrix a2   (first derivative)
a2 = -2*Diagonal(y) * d1;

# assemble matrix a3   (constant term)
y2c_inv = 1.0./y2c;
a3 = -m^2*Diagonal(y2c_inv);

# now combine into A0, A1 and A2
A0 = convert(Array{Complex{Float64},2},a1+a2+a3);
A1 = convert(Array{Complex{Float64},2},m*Diagonal(ones(n))/M);
A2 = convert(Array{Complex{Float64},2},Diagonal(y2)/M);


return A0,A1,A2

end

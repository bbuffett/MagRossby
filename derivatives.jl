function derivatives(y)

# evaluates the matrices for the first and second derivatives
# with respect to y
#
#
# input
#  y -  meridional coordinate

#
# output
# d1, d2 - matrices for first and second derivative

dy = y[2]-y[1];
ngrid = length(y);
a = Vector{Float64}(undef,ngrid);
b = Vector{Float64}(undef,ngrid-1);
c = Vector{Float64}(undef,ngrid-1);

#
# second derivative
a[1:ngrid] .= -2.0/dy^2;
b[1:ngrid-1] .= 1.0/dy^2;
d2 = SymTridiagonal(a,b);

# first derivative
c[1:ngrid-1] .= 1.0/(2*dy);
d1 = Tridiagonal(-c,zeros(ngrid),c);


return d1,d2

end

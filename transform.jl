function transform(y,bxp,byp)

# This function transforms the modified variable
# to the original components of the magnetic
# perturbation according to Eq (13) of B&M (2019)
# We also impose the conditions at the pole and
# interpolate the solution to give better resolution
# in latitude

# input
#  y        - grid points in meridional coordinate
#  bxp,byp  - modified magnetic perturbation
#
# output
# ynew   - updated grid points including pole at y=1
# bx,by  - original magnetic perturbations


# apply transform back to original variables
y2 = y.*y;
yc2 = -(y2 .-1.0);
f = sqrt.(yc2);
bx = bxp .* f;
by = byp ./ f;

# impose boundary condition at pole y = 1
last_lat = asin(y[end]);
dlat = (pi/2 - last_lat)/4;
lat_intpol = collect(last_lat : dlat : pi/2);
y_intpol = sin.(lat_intpol);
yc = -(y_intpol.*y_intpol .- 1.0);
yc_end = 1.0 - y[end]^2;
by_intpol = by[end] * (yc / yc_end);
bx_intpol = bx[end] * (yc / yc_end);

# append extrapolation to solution
ynew = copy(y);
append!(ynew,y_intpol[2:end]);
append!(bx,bx_intpol[2:end]);
append!(by,by_intpol[2:end]);


return ynew,bx,by

end

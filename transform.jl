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

# append grid point at the pole y=1
 ynew = copy(y);
 push!(ynew,1.0);

# append solution at the pole
  if (m > 1)
    push!(bx,0.0);
    push!(by,0.0);
  else
    bxend = 2*bx[end]-bx[end-1];  # interpolate
    byend = 2*by[end]-by[end-1];
    push!(bx,bxend);
    push!(by,byend);
  end
    


return ynew,bx,by

end

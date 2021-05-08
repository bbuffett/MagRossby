function mgrid(ngrid)

# evaluates the grid in the meridional coordinate y=cos(theta)
# and also returns the corresponding values of latitude. The
# end points at the equator and the pole are excluded.
#
#
# input
# ngrid  - number of grid points
#
# output
# y  - grid values
# lat - grid values in latitude (radian)
#


# Evaluates the grid
dy = 1.0 / (ngrid+1);
y = collect(dy : dy : 1.0-dy);

# corresponding latitude (degree)
lat = asin.(y)*180/pi;

return y, lat

end

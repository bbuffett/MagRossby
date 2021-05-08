function divergence(t,colat,long,vph,vth,aph,ath)

# computes the divergence of the velocity and acceleration
# on a lat,long grid at a series of snapshots.
#
# input
#
# t -  array of times
# colat  - array of colatitudes
# long  - array of longitudes
# vph,vth  - phi and theta components of velocity (km/yr)
# aph,ath  - phi and theta components of acceleration (km/yr^2)

#
# output
#
# vdiv   -  divergence of velocity
# adiv   -  divergence of acceleration


# allocate memory for divergence
vdiv = Array{Float64,3}(undef,length(t),length(colat),length(long));
adiv = Array{Float64,3}(undef,length(t),length(colat),length(long));

# radius of core-mantle boundary
rcmb = 3480.0;  #  in km

# loop over times
for i = 1 : length(t)

#  loop over latitude and longitude
   nlat = length(colat)
   nlong = length(long)
   for j = 1 : nlat

      rsinth = rcmb * sin(colat[j]);
      if (j < 0.1*nlat) || (j > nlat - 0.1*nlat)
         f = sin(colat[j]);
      else
         f = 1.0;
      end

      for k = 1 : nlong

         # derivative of vth wrt colatitude
         if (j == 1)
           dlat = colat[j+1];
           dvth = f*(sin(colat[j+1])*vth[i,j+1,k] -
                     0.0) / (rsinth * dlat);
           dath = f*(sin(colat[j+1])*ath[i,j+1,k] -
                     0.0) / (rsinth * dlat);

         elseif (j == nlat )
            dlat = pi - colat[j-1];
            dvth = f*(0.0 - sin(colat[j-1])*vth[i,j-1,k] ) /
                    (rsinth * dlat);
            dath = f*(0.0 - sin(colat[j-1])*ath[i,j-1,k] ) /
                    (rsinth * dlat);

         else
            dlat = colat[j+1] - colat[j-1];
            dvth = f*(sin(colat[j+1])*vth[i,j+1,k] -
                  sin(colat[j-1])*vth[i,j-1,k])/(rsinth*dlat);
            dath = f*(sin(colat[j+1])*ath[i,j+1,k] -
                  sin(colat[j-1])*ath[i,j-1,k])/(rsinth*dlat);
         end

         # derivative of vphhi wrt longitude
        if (k == 1)
            dlong = long[k+1] - (long[nlong]-2*pi);
            dvph = f*(vph[i,j,k+1] - vph[i,j,nlong])/(rsinth * dlong);
            daph = f*(aph[i,j,k+1] - aph[i,j,nlong])/(rsinth * dlong);

         elseif (k == nlong)
            dlong = long[1] + 2*pi - long[k-1];
            dvph = f*(vph[i,j,1] - vph[i,j,k-1])/(rsinth * dlong);
            daph = f*(aph[i,j,1] - aph[i,j,k-1])/(rsinth * dlong);

         else
            dlong = long[k+1] - long[k-1];
            dvph = f*(vph[i,j,k+1] - vph[i,j,k-1])/(rsinth * dlong);
            daph = f*(aph[i,j,k+1] - aph[i,j,k-1])/(rsinth * dlong);
         end

         # divergence
         vdiv[i,j,k] = dvth + dvph;
         adiv[i,j,k] = dath + daph;

      end

   end

end

return vdiv, adiv

end

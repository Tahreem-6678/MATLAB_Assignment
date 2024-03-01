function distance = sat_dist(h,e)
%This programme calculates the distance to a satellite.
%
%  Inputs:
%       h = altitude in km above the earth's surface
%       e = elevation in degrees, as seen by the observer on earth.

%  output:
%       distance in km from the observer to the satellite.

re = 6370; %Earth radius in km.
z=re.*sind(e);
distance = sqrt(z.^2+h.^2+2*re.*h)-z;
end


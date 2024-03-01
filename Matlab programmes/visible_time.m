function tvis = visible_time(h,e)
% This function assumes an object flies overhead in a circular orbit ,,
% passing from e degrees above one horizon to e degrees above the diametrically opposite horizon 
% via the zeniith.  It calculates the time in seconds for which the object is visible.
% The earth's rotation is neglected.
%   Inputs
%   h   Height of the object in km.
%   e   Minimum elevation in degrees
G= 6.674e-11; %Universal gravitational const.
M=5.972e24; %Mass of earth.
er=6371;

period=2*pi*sqrt(((h+6370)^3*1e9)/(G*M)) %Orbital period in seconds.
d = sat_dist(h,e);  %Distance between satellite and observer at lowest visible elev.
ang=acos((2*er*(er+h)+h.^2-d.^2)/(2*er*(er+h))); %The angle subtended at earth's centre by 
% the line joining sat at its lowerst visible elevation to earth station.
tvis=ang*period/pi
end
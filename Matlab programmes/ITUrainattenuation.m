function Ap = ITUrainattenuation(frequency, theta, h0, hS, R001, phi, percentage);

% Function to calculate the long-term rain attenuation
% statistics from point rainfall rate using ITU-R P618.7

% sample run to plot attenuation against frequency in range 1-55 GHz:
% plot([1:55],ITUrainattenuation([1:55], 45, 4.077, 0.64, 43, 37.229, 0.01))

% Verification of Pratt and Bostian example (page 338-339):
% ITUrainattenuation(28, 45, 4.077, 0.64, 43, 37.229, 0.01)

% Note: No more than one input value can be a vector.

% Inputs:
%  frequency: carrier frequency (GHz)
%  theta: elevation angle (degrees)
%  h0: mean zero degree isotherm height above mean sea-level (km)
%  hS: height above mean sea level of the earth station (km)
%  R001: point rainfall rate for a location for 0.01% of an average year (mm/h).
%        This can be obtained from ITU-R P.837 maps
%  phi: latitude of the earth station (degrees)
%  percentage: percentage of an average year to convert atenuation to (0.01 is default)

% Output:
%  Ap: Estimated Rain Attenuation (dB)

% Intermediate values calculated by the function:
%  hR: the mean rain height above mean sea-level (km) 
%  hR = h0 + 0.36, 
%       where h0 is the mean zero degree isotherm height above mean sea-level
%       Based on ITU-R P.839, 1<h0<4.5 km
%  LE: effective path length (km)
%  Ls: slant-path length, below the rain height (km)
%  LG: the horizontal projection of the slant-path length (km)

% Check the input parameters for sensible values
if length(find(R001>120 | R001 < 1))~=0
    fprintf('ERROR - the rain rate is out of bounds\n');
    return
end

if length(find(abs(phi)>90))~=0
    fprintf('ERROR - the earth station latitude range is out of bounds\n');
    return
end

if length(find(frequency>55 | frequency < 0.003))~=0
    fprintf('ERROR - the frequency range is out of bounds\n');
    return
end

if length(find(theta>90 | theta < 0))~=0
    fprintf('ERROR - the elevation angle is out of bounds\n');
    return
end

Re = 8500;  % effective radius of the earth(km)
theta = theta*(pi/180);	% convert theta from degrees to radians

%  *** Step-1 ***
% Compute rain height, hR
hR = h0 + 0.36; % in km

% check that the station height is not above the rain height)
if length(find(hS-hR>0))~=0
	fprintf('ERROR - station height above rain height\n');
	return
end

%  *** Step-2 ***
% Compute slant path length, Ls
if  (theta >= 5*pi/180)	% converted 5 degrees to radians for comparison with theta
    Ls = (hR-hS)./sin(theta);
else  
    nr1 = 2*(hR-hS);
    dr1 = (sqrt((sin(theta).^2)+((2*(hR-hS))./Re))) + sin(theta);
    Ls = nr1./dr1;
end

%  *** Step-3 ***
% Calculate the horizontal projection, LG, of the slant-path length
LG = Ls .* cos(theta);

%  *** Step-4 ***
% Obtain the rainfall rate, R001, exceeded for 0.01% of an average year
% - obtained from ITU-R P.837 maps (in mm/hr)

%  *** Step-5 ***
% Obtain the specific attenuation, gammaR, using ITU-R P.838-2
% This has been determined only for cases of circular polarization, where
% the polarization tilt angle relative to the horizontal (tau) = 45 degrees
% => the cos(2*tau) terms in equations 4 & 5 = 0

% ****************************************************************************
% Calculating the values of kH, kV, alphaH and alphaV: see ITU-R P.838-2
% H->horizontal polarization
% V->vertical polarization

% calculating kH
% coefficients for horizontal polarization
aHj=[0.3364,0.7520,-0.9466];
bHj=[1.1274,1.6644,2.8496]; 
cHj=[0.2916,0.5175,0.4315];
mkH=1.9925;
ckH=-4.4123;

klogH=0;
for j=1:3
       klogH = klogH + aHj(j).*exp(-((log10(frequency)-bHj(j))./cHj(j)).^2);
end
klogH = klogH + mkH.*log10(frequency) + ckH;
kH=10.^(klogH);

%calculating kV
% coefficients for vertical polarization
aVj=[0.3023,0.7790,-1.0022];
bVj=[1.1402,1.6723,2.9400]; 
cVj=[0.2826,0.5694,0.4823];
mkV=1.9710;
ckV=-4.4535;

klogV=0;
for j=1:3
       klogV = klogV + aVj(j).*exp(-((log10(frequency)-bVj(j))./cVj(j)).^2);
end
klogV = klogV + mkV.*log10(frequency) + ckV;
kV=10.^(klogV);

% calculating alphaH
% coefficients for horizontal polarization
aHi=[0.5564,0.2237,-0.1961,-0.02219];
bHi=[0.7741,1.4023,0.5769,2.2959];   
cHi=[0.4011,0.3475,0.2372,0.2801];
malphaH=-0.0816;
calphaH=0.8993;

alphaH=0;
for i=1:4
       alphaH = alphaH + aHi(i).*exp(-((log10(frequency)-bHi(i))./cHi(i)).^2);
end
alphaH = alphaH + malphaH.*log10(frequency) + calphaH;

% calculating alphaV 
% coefficients for vertical polarization
aVi=[0.5463,0.2158,-0.1693,-0.01895];
bVi=[0.8017,1.4080,0.6353,2.3105]; 
cVi=[0.3657,0.3636,0.2155,0.2938];
malphaV=-0.07059;
calphaV=0.8756;

alphaV=0;
for i=1:4
       alphaV = alphaV + aHi(i).*exp(-((log10(frequency)-bVi(i))./cVi(i)).^2);
end
alphaV = alphaV + malphaV.*log10(frequency) + calphaV;

k = (kH + kV)./2;
alpha = (kH.*alphaH + kV.*alphaV)/(2.*k); 

gammaR = k.*((R001).^alpha);

% ****************************************************************************

%  *** Step-6 *** 
% Calculate the horizontal reduction factor, r001, for 0.01% of the time
dr2 = (0.78.*sqrt(LG.*gammaR./frequency)) - (0.38.*(1-exp((-2).*LG)));
r001 = 1./(1 + dr2);

%  *** Step-7 ***
% Calculate the vertical adjustment factor, v001, for 0.01% of the time
zeta = atan((hR - hS)./(LG.*r001)); % zeta is in radians

if zeta > theta % note - theta is in radians
    LR = (LG.*r001)./(cos(theta)); % in km
else 
    LR = (hR - hS)./(sin(theta)); % in km  
end

%if abs(phi) < 36	% note - phi is in degrees
%    chi = 36 - abs(phi);
%else
%    chi = 0; % in degrees
%end
temp=find(abs(phi)<36);
chi(temp)=36-abs(phi(temp));
temp=find(abs(phi)>=36);
chi(temp)=0;

dr3 = sqrt(sin(theta)).*(31*(1-exp(-((theta*180/pi)./(1+chi)))).*((sqrt(LR.*gammaR))./frequency.^2)-0.45);
v001 = 1./(1+dr3);
  
%  *** Step-8 ***
% Calculate the effective path length, LE 
LE = LR.*v001; % in km

%  *** Step-9 ***
% Obtain the predicted attenuation exceeded for 0.01% of an average year 
A001 = gammaR.*LE; % in dB
  
% Finally calculating the estimated attenuation to be exceeded for other percentages of an
% average year (in the range 0.001% to 5%)

if percentage>=1 | abs(phi)>=36
   beta=0;
elseif percentage<1 & abs(phi)<36 & theta>=25*180/pi
   beta=-0.005.*(abs(phi)-36);
else 
   beta=(-0.005.*(abs(phi)-36)) + 1.8 - (4.25.*sin(theta));
end

exp_term = -(0.655+(0.033.*log(percentage))-(0.045.*log(A001))-(beta.*(1-percentage).*sin(theta)));
Ap = A001.*((percentage./0.01).^exp_term);

% end of program




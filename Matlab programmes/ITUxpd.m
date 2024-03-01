function XPDp = ITUxpd(frequency, tau, theta, Ap, percentage)

% Function to calculate the loss due to the long-term statistics of
% hydrometeor-induced cross-polarization using ITU P.618-7
%
% Inputs:
%  frequency: frequency (GHz)
%  tau: tilt angle of the linearly polarized electric field vector with
%       respect to the horizontal (tau = 45 degrees for circular polarization)
%  theta: path elevation angle (degrees)
%  Ap: rain attenuation (dB) exceeded for the required percentage of time, p
%  percentage: percentage of an average year (1, 0.1, 0.01, or 0.001%)

% Note: Any one of the inputs, except percentage, can be a vector.

% Outputs:
%  XPDp: cross-polarization discrimination (dB)

if length(find(frequency > 35 | frequency < 8))~=0
   fprintf('ERROR: The frequency range is out of bounds\n');
   return;
end

if length(find(theta > 60))~=0
   fprintf('ERROR: The path elevation angle is out of bounds\n');
   return;
end

if percentage~=1 & percentage~=0.1 & percentage~=0.01 & percentage~=0.001
	fprintf('ERROR - percentage must be 1, 0.1, 0.01 or 0.001\n')
	return
end

% convert angles to radians 
theta = theta.*(pi/180);
tau = tau.*(pi/180);

%  *** Step-1 ***
% Calculating the frequency-dependent term:
Cf = 30.*log10(frequency);
   
%  *** Step-2 ***
% Calculating the rain attenuation dependent term:
temp=find(frequency >= 8 & frequency <= 20);
if length(temp)~=0
   Vf(temp) = 12.8.*(frequency(temp).^(0.19));
end
  
temp=find(frequency > 20 & frequency <= 35);
if length(temp)~=0
   Vf(temp) = 22.6;
end
  
CA = Vf.*log10(Ap);
  
%  *** Step-3 ***
% Calculating the polarization improvement factor:
Ctau = -10.*log10(1 - 0.484.*(1 + cos(4.*tau)));
 
%  *** Step-4 ***
% Calculating the elevation angle-dependent term:
Ctheta = -40.*log10(cos(theta));
 
%  *** Step-5 ***
% Calculating the canting angle dependent term: sigma = effective
% standard deviation of the raindrop canting angle distribution (degrees)
if percentage == 1
   sigma = 0;
elseif percentage == 0.1
   sigma = 5;
elseif percentage == 0.01
   sigma = 10;
elseif percentage == 0.001
   sigma = 15;
end

Csigma = 0.0052.*(sigma.^2);

%  *** Step-6 *** 
% Calculating rain XPD not exceeded for p% of the time:
XPDrain = Cf - CA + Ctau + Ctheta + Csigma; % (dB)
 
%  *** Step-7 *** 
% Calculating the ice crystal dependent term:
Cice = XPDrain.*(0.3 + 0.1.*log10(percentage))./2;

%  *** Step-8 ***   
% Calculating the XPD not exceeded for p% of the time, including the effects of ice:
XPDp = XPDrain - Cice;
  
% end of program

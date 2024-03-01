function Gain = ITUdiversitygain(frequency,d,A,theta,psi)

% Function to calculate the Site Diversity Gain (dB) using ITU-R P.618-7

% Sample run:
% plot([10:30],ITUdiversitygain([10:30],5,25,45,90))

% Inputs:    
%  frequency: frequency (GHz): frequency range of applicability: 10 -30 GHz
%  d: separation between the two sites (km)
%  A: Single site rain attenuation (dB). Could use ITU model.
%  theta: path elevation angle (degrees)
%  psi: angle (degrees) made by the azimuth of the propagation path with 
%       respect to the baseline between sites, chosen such that psi<=90

% Output:
%  Gain: net diversity gain (dB)

if length(find(abs(psi)>90 | abs(psi)<0))~=0
    fprintf('ERROR - the baseline angle, psi, is out of bounds\n');
    return
end

if length(find(frequency>30 | frequency < 10))~=0 
    fprintf('ERROR - the frequency range is out of bounds\n');
    return
end

if length(find(A<0))~=0
	fprintf('ERROR - the rain attenuation must be positive\n')
	return
end

% *** Step-1 ***	
% Calculate the gain contributed by the spatial separation:

a = (0.78.*A) - (1.94.*(1-exp(-0.11.*A)));
b = 0.59.*(1-exp(-0.1.*A));

Gd=a.*(1-exp(-b.*d));

% *** Step-2 ***
% Calculate the frequency-dependent gain:
Gf = exp(-0.025.*frequency);

% *** Step 3 ***	
% Calculate the gain term dependent on elevation angle:
Gtheta = 1 + (0.006.*theta);

% *** Step-4 ***
% Calculate the baseline-dependent term:
Gpsi = 1 + (0.002.*psi);

% *** Step-5 ***
% Compute the net diversity gain as the product:
Gain = Gd.*Gf.*Gtheta.*Gpsi;	% (dB)

% end of program


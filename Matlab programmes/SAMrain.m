function absorption=SAMrain(freqGHz, elevation, rainrate, Hstation, Hice)

% function to calculate the absorption due to hydrometeors according to
% the SAM model (Pratt and Bostian)
%	Written by AJS (11/11/04)
%
% inputs:
%	freqGHz : The wave frequency in GHz
%	elevation : Elevation angle to the satellite (degrees)
%	rainrate : The rain rate (mm/h) [scalar]
%	Hstation : The height of the ground station (km) [must be below Hstorm]
%	Hice : The height of the zero degree isotherm (km), if set to a
%		negative number this is the magnitude of the latitude.
%
% Note: any ONE of the input values can be a vector except rainrate which must
% be a scalar.
%
% output
%	absorption : absorption in dB
%
% sample call:
%	plot([10:50],SAMrain([10:50], 45, 43, 0.64, -37.229))
%
% Will plot absorption as a function of frequency for Blacksburg, Virginia
% (see example in Pratt and Bostian)

% check input values
	if length(find(freqGHz<8.5 | freqGHz>164))~=0
		fprintf('ERROR - frequency must be between 8.5 and 164 GHz\n')
		return
	end

	if length(find(elevation<0 | elevation>90))~=0
		fprintf('ERROR - elevation must be between 0 and 90 degrees\n')
		return
	end

	if rainrate<0
		fprintf('ERROR - rainrate must be positive\n')
		return
	end
		
% convert elevation to radians
	elevation=elevation*pi/180;
	
% find constants of specific attenuation
	temp=find(freqGHz>=2.9 | freqGHz<=54);
	if length(temp)~=0
		a(temp)=4.21e-5*freqGHz(temp).^2.42;
	end
	temp=find(freqGHz>54 & freqGHz<180);
	if length(temp)~=0
		a(temp)=4.09e-2*freqGHz(temp).^0.699;
	end
	temp=find(freqGHz>=8.5 | freqGHz<=25);
	if length(temp)~=0
		b(temp)=1.41*freqGHz(temp).^-0.0779;
	end
	temp=find(freqGHz>25 & freqGHz<164);
	if length(temp)~=0
		b(temp)=2.63*freqGHz(temp).^-0.272;
	end

% find ice height and then storm height
% for the Hice values, if they are negative this means that it
% is a latitude and the SAM formulas should be used.
	temp=find(Hice<0 & Hice>=-30);
	if length(temp)~=0
		Hice(temp)=4.8;
	end
	temp=find(Hice<-30 & Hice>=-68);
	if length(temp)~=0
		Hice(temp)=7.8-0.1*abs(Hice(temp));
	end
	temp=find(Hice<-68);
	if length(temp)~=0
		Hice(temp)=1.0;	% this is a high-latitude correction to the model
	end

	if rainrate<=10;
		Hstorm=Hice;
	elseif rainrate>10
		Hstorm=Hice+log10(rainrate/10);
	end

% check that station height is not above storm height
	if length(find(Hstation>Hstorm))~=0
		fprintf('ERROR - Station height is above storm height\n')
		return
	end
	
% calculate path length in rain
	L=(Hstorm-Hstation)./sin(elevation);
	
% calculate absorption
	gamma=1/22;	% empirical constant
	if rainrate<=10
		absorption=a.*L.*rainrate.^b;
	elseif rainrate>10
		absorption=a.*rainrate.^b.*(1-exp(-gamma.*b*log(rainrate/10).*L ...
			.*cos(elevation)))./(gamma*b*log(rainrate/10).*cos(elevation));
	end
	

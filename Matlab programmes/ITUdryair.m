function gamma0 = ITUdryair(frequency,pressure,temperature);

% Function to calculate the specific attenuation due to dry air in the
% frequency range 1 to 350 GHz (ITU-R P.676-5, Annex 2)

% Sample run:
% loglog([1:350],ITUdryair([1:350],1013,15)) reproduces the curve
% on page 15 in ITU-R P676-5

% Inputs:
%  frequency: carrier frequency (GHz)
%  pressure: pressure (hPa)
%  temperature: temperature (degrees C)
%             - annual mean surface temperature values can be obtained from maps given in
%               ITU-R P.1510, when no adequate temperature data is available

% Output:
%  gamma0: Specific attenuation due to dry air (dB/km)

% Intermediate values calculated by the function:
%  rp: pressure-dependent coefficient = pressure/1013; 
%  rt: temperature-dependent coefficient = 288/(273+temperature);
%  Coefficients: gamma054, gamma054a, gamma057, gamma060, gamma063, gamma066, gamma066a
%              - eeta1, eeta2, a, b, c, d, zeta1, zeta2

if length(find(frequency<1 | frequency>350))~=0
   fprintf('Error: frequency is out of bounds.\n');
   return
end

if length(find(temperature<-10 | temperature>25))~=0
   fprintf('Error: Surface temperature is out of bounds.\n');
   return
end

rp = pressure./1013;
rt = 288./(273 + temperature);

gamma054 = 2.136.*(rp.^1.4975).*(rt.^(-1.5852)).*(exp(-2.5196.*(1 - rt)));
gamma054a = 2.128.*(rp.^1.4954).*(rt.^(-1.6032)).*(exp(-2.5280.*(1 - rt)));
gamma057 = 9.984.*(rp.^0.9313).*(rt.^2.6732).*(exp(0.8563.*(1 - rt)));
gamma060 = 15.42.*(rp.^0.8595).*(rt.^3.6178).*(exp(1.1521.*(1 - rt)));
gamma063 = 10.63.*(rp.^0.9298).*(rt.^2.3284).*(exp(0.6287.*(1 - rt)));
gamma066 = 1.944.*(rp.^1.6673).*(rt.^(-3.3583)).*(exp(-4.1612.*(1 - rt)));
gamma066a = 1.935.*(rp.^1.6657).*(rt.^(-3.3714)).*(exp(-4.1643.*(1 - rt)));

eeta1 = (6.7665.*(rp.^(-0.5050)).*(rt.^0.5106).*(exp(1.5663.*(1 - rt)))) - 1;
eeta2 = (27.8843.*(rp.^(-0.4908)).*(rt.^0.8491).*(exp(0.5496.*(1 - rt)))) - 1;
zeta1 = (6.9575.*(rp.^(-0.3461)).*(rt.^0.2535).*(exp(1.3766.*(1 - rt)))) - 1;
zeta2 = (42.1309.*(rp.^(-0.3068)).*(rt.^1.2023).*(exp(2.5147.*(1 - rt)))) - 1;

a = log(eeta2./eeta1)./log(3.5);
b = (4.^a)./eeta1;
c = log(zeta2./zeta1)./log(3.5);
d = (4.^c)./zeta1;

if length(find(frequency<=54))~=0
   f1 = frequency(find(frequency<=54));
   gamma0_a = (((7.34.*(rp.^2).*(rt.^3))./((f1.^2)+(0.36.*(rp.^2).*(rt.^2)))) + ((0.3429.*b.*gamma054a)./...
              (((54-f1).^a)+b))).*(f1.^2).*(10.^(-3));
else
   gamma0_a=[];          
end

if length(find(frequency>54 & frequency<=66))~=0
    if length(find(frequency>54 & frequency<=60))~=0
       N = 0; 
       f2a = frequency(find(frequency>54 & frequency<=60));
       gamma0_b1 = exp((((54.^(-N).*log(gamma054).*(f2a-57).*(f2a-60).*(f2a-63).*(f2a-66))./1944)-((57.^(-N).*log(gamma057).*(f2a-54).*(f2a-60).*(f2a-63).*(f2a-66))./486)...
                   +((60.^(-N).*log(gamma060).*(f2a-54).*(f2a-57).*(f2a-63).*(f2a-66))./324)-((63.^(-N).*log(gamma063).*(f2a-54).*(f2a-57).*(f2a-60).*(f2a-66))./486)...
                   +((66.^(-N).*log(gamma066).*(f2a-54).*(f2a-57).*(f2a-60).*(f2a-63))./1944)).*(f2a.^N));
    else    
       gamma0_b1 = [];               
	end
	
	if length(find(frequency>60 & frequency<66))~=0
       N = -15; 
       f2b = frequency(find(frequency>60 & frequency<66));
       gamma0_b2 = exp((((54.^(-N).*log(gamma054).*(f2b-57).*(f2b-60).*(f2b-63).*(f2b-66))./1944)-((57.^(-N).*log(gamma057).*(f2b-54).*(f2b-60).*(f2b-63).*(f2b-66))./486)...
                   +((60.^(-N).*log(gamma060).*(f2b-54).*(f2b-57).*(f2b-63).*(f2b-66))./324)-((63.^(-N).*log(gamma063).*(f2b-54).*(f2b-57).*(f2b-60).*(f2b-66))./486)...
                   +((66.^(-N).*log(gamma066).*(f2b-54).*(f2b-57).*(f2b-60).*(f2b-63))./1944)).*(f2b.^N));
    else    
       gamma0_b2 = [];          
    end
     gamma0_b = [gamma0_b1,gamma0_b2];
else    
    gamma0_b = [];    
end;    

if length(find(frequency>=66 & frequency<120))~=0
   f3 = frequency(find(frequency>=66 & frequency<120));
   gamma0_c = ((0.2296.*d.*gamma066a)./(((f3-66).^c) + d) + ((0.286.*(rp.^2).*(rt.^3.8))./(((f3-118.75).^2) ...
              + (2.97.*(rp.^2).*(rt.^1.6))))).*(f3.^2).*(10.^(-3));
else    
   gamma0_c = [];        
end
    
if length(find(frequency>=120 & frequency<=350))~=0
   f4 = frequency(find(frequency>=120 & frequency<=350));
   gamma0_d = ((3.02.*(10.^(-4)).*(rp.^2).*(rt.^3.5)) + ((1.5827.*(rp.^2).*(rt.^3))./(f4-66).^2) + ...
              ((0.286.*(rp.^2).*(rt.^3.8))./(((f4-118.75).^2) + (2.97.*(rp.^2).*(rt.^1.6))))).*(f4.^2).*(10.^(-3));           
else    
    gamma0_d = [];          
end        

gamma0=[gamma0_a,gamma0_b,gamma0_c,gamma0_d];

% end of program

function loss=freespace(frequency, distance)
% function to calculate the free space loss (following Shankar, pg13-14)
%
% Inputs:
%	frequency :	Carrier frequency (in MHz)
%	distance :	distance between tx and rx (in km, should be >1)
%
% Ouputs :
%	loss : Signal loss (in dB)

% warn if the any distances are below 1 km (the model will still give an
% output, you will have been warned that it may be in error)
if length(find(distance<1))~=0
	fprintf('WARNING - distance set below 1km in freespace loss\n');
end

loss=32.44 + 20*log10(frequency) + 20*log10(distance);


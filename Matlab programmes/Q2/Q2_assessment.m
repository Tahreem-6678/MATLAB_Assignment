
% Satellite parameters
h = 400; % Altitude (km)
e_min = 10; % Minimum elevation (degrees)

% Earth radius
er = 6371; % km

% Data rate
data_rate = 10 * 1e3; % kbit/s (convert to kb/s)

% Calculate distance at minimum elevation
distance = sat_dist(h, e_min);

% Calculate free space loss
loss = freespace(5e9, distance); % Frequency in GHz, distance in km

% Display free space loss
fprintf("Free space loss when data collection start: %.2f dB\n", loss);

% Calculate orbital period
G = 6.674e-11; % Gravitational constant
M = 5.972e24; % Earth mass
period = 2*pi*sqrt(((h + er)^3 * 1e9) / (G * M)); % Orbital period in seconds

% Calculate visible time
tvis = visible_time(h, e_min);

% Calculate minimum required data rate
required_rate = data_rate / (tvis / 1e3); % Convert to kb/s

% Display required data rate
fprintf("Required data rate to download data: %.2f kb/s\n", required_rate);

% Simplified function for oxygen attenuation (adapted from ITUdryair.m)
function att_oxygen = oxygen_att(frequency, p_ground, T_ground, h)
    % Reference pressure and temperature
    p_ref = 1013; % hPa
    T_ref = 273 + 15; % K (assuming 15°C at 10 km)

    % Scale height for oxygen (assumed constant)
    H = 8.5; % km

    % Oxygen attenuation coefficient at reference conditions
    alpha_ref = 0.12; % dB/km (assumed value)

    % Pressure at satellite altitude
    p_sat = p_ref * exp(-h / H);

    % Temperature adjustment factor
    T_factor = T_ref / (273 + T_ground);

    % Oxygen attenuation
    att_oxygen = alpha_ref * p_sat / p_ref * T_factor * (frequency / 50e9)^2;

end

% **Satellite parameters**
h = 500; % km
eirp = 100; % W
frequency = 50 * 1e9; % Hz (convert GHz to Hz)
dish_diameter = 0.9; % m

% **Ground station parameters**
p_ground = 1013; % hPa (assumed)
T_ground = 20; % °C

% **Minimum required Eb/No**
Eb_No_min = -120 - 30; % dBW/Hz (assuming 30 dB noise figure)

% **Boltzmann's constant**
kB = 1.38e-23; % J/K

% **a. Scale height for oxygen pressure**
fprintf("Scale height for oxygen pressure: %.2f km\n", H);

% **b. Oxygen attenuation when overhead (e=90°)**
att_overhead = oxygen_att(frequency, p_ground, T_ground, h);
fprintf("Oxygen attenuation (overhead): %.2f dB/km\n", att_overhead);

% **c. Oxygen loss vs. elevation angle**
elev_angles = 10:1:90; % degrees
loss_oxygen = att_overhead ./ sind(elev_angles);

% Plot oxygen loss vs. elevation angle
figure(1);
plot(elev_angles, loss_oxygen);
xlabel('Elevation Angle (degrees)');
ylabel('Oxygen Loss (dB)');
title('Oxygen Loss vs. Elevation Angle');
grid on;

% **d. Receive antenna gain**
A_r = 10 * log10(pi * dish_diameter^2 / (4 * pi));
fprintf("Receive antenna gain: %.2f dBi\n", A_r);

% **e. Distance and free-space loss vs. elevation angle**
distances = sat_dist(h, elev_angles);
loss_free_space = freespace(frequency, distances);

% Plot distance vs. elevation angle
figure(2);
plot(elev_angles, distances);
xlabel('Elevation Angle (degrees)');
ylabel('Distance (km)');
title('Distance vs. Elevation Angle');
grid on;

% Plot free-space loss vs. elevation angle
figure(3);
plot(elev_angles, loss_free_space);
xlabel('Elevation Angle (degrees)');
ylabel('Free-Space Loss (dB)');
title('Free-Space Loss vs. Elevation Angle');
grid on;

% **f. Total loss and minimum elevation angle**
loss_total = loss_oxygen + loss_free_space;

% Plot total loss vs. elevation angle
figure(4);
plot(elev_angles, loss_total);
xlabel('Elevation Angle (degrees)');
ylabel('Total Loss (dB)');
title('Total Loss vs. Elevation Angle');
grid on;

% Calculate maximum allowed loss based on link budget
max_loss = Eb_No_min + 10 * log10(kB * T_ground / (eirp * dish_diameter^2)); % dB

% **g. Minimum elevation angle for reception**
min_elev_angle = find(loss_total <= max_loss, 1, 'first');

fprintf("Minimum elevation angle for reception: %d degrees\n", min_elev_angle);

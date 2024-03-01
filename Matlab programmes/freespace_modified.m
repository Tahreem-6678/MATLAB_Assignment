% using function loss=freespace(frequency, distance) storeed in freespace.m

% Generate data for plotting
frequencies = [1 2]; % GHz
distances = 100:1:1000; % km

% Calculate loss for each frequency and distance
loss_data = zeros(length(frequencies), length(distances));
for i = 1:length(frequencies)
    loss_data(i,:) = freespace(frequencies(i) * 1e3, distances); % Convert GHz to MHz
end

% Plot the results
figure(1); % Set figure number
hold on;
for i = 1:length(frequencies)
    plot(distances, loss_data(i,:), 'DisplayName', [num2str(frequencies(i)) ' GHz']);
end

% Customize the plot
grid on;
legend;
xlabel('Distance (km)');
ylabel('Free Space Loss (dB)');
title('Free Space Loss vs. Distance for Different Frequencies');

% Add text description
text(200, 140, 'As distance increases, free space loss increases.', 'FontSize', 12, 'Color', 'red');

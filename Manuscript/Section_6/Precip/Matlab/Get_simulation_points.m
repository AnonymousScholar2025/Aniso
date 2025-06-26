%mesh
tv = load('delTV.txt');
loc = load('delLoc.txt');
border = load('delBorder.txt');
% Remove NaN from border
border2 = border(~any(isnan(border), 2), :);

% Sampling from box [a,b]x[c,d]
a = max(min(loc(:, 1)), min(border2(:, 1)));
b = min(max(loc(:, 1)), max(border2(:, 1)));
c = max(min(loc(:, 2)), min(border2(:, 2)));
d = 7250;
rng(123);
n_samples = 13000; % Number of points in mesh to simulate from
observation_locations = [a + (b-a).*rand(n_samples, 1), c + (d-c).*rand(n_samples, 1)];

% Select points within border
mask = inpolygon(observation_locations(:, 1), observation_locations(:, 2), border(:, 1), border(:, 2));
obs_locations = observation_locations(mask, :);

% Save locations and height
writematrix(obs_locations, 'sim_data/observation_locations_sim.txt');

% Load simulated points
obs_locations = readmatrix('sim_data/observation_locations_sim.txt');


% Plot the border
plot(border(:,1), border(:,2));
hold on;

% Add the observation locations
scatter(obs_locations(:, 1), obs_locations(:, 2));

hold off;
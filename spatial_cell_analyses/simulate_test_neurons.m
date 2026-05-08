% =========================================================================
% SYNTHETIC SPATIAL CELL GENERATOR
% Simulates a Trajectory, Place Cell, Grid Cell, and Head Direction Cell
% =========================================================================

%% 1. Setup Parameters
% Time and Space
fs = 50;                  % Sampling rate (Hz)
dt = 1/fs;                % Time step
duration = 600;           % Total time (seconds)
t = 0:dt:duration;        % Time vector
N = length(t);
boxSize = 100;            % Size of the square arena (cm)

% Place Cell Parameters
pc_centers = [25, 75; 80, 20]; % Two fields at (X,Y)
pc_width = 10;                 % Gaussian standard deviation
pc_max_rate = 20;              % Peak firing rate (Hz)

% Grid Cell Parameters
grid_spacing = 30;             % Distance between grid nodes (cm)
grid_orientation = pi/6;       % Tilt of the grid
grid_max_rate = 15;            % Peak firing rate (Hz)

% Head Direction (HD) Cell Parameters
hd_pref = pi/2;                % Preferred direction (90 degrees / North)
hd_kappa = 3;                  % Tuning width (higher = narrower)
hd_max_rate = 30;              % Peak firing rate (Hz)

%% 2. Generate Animal Trajectory (Correlated Random Walk)
x = zeros(1, N);
y = zeros(1, N);
hd = zeros(1, N); % Head direction

% Start in the center of the box
x(1) = boxSize / 2;
y(1) = boxSize / 2;
hd(1) = rand() * 2 * pi; % Random starting angle

% Locomotion Parameters
speed = 20;       % Constant forward speed (cm/s)
turn_std = 0.2;   % How much the animal turns per time step (radians)
                  % (Lower = straighter paths, Higher = tighter circles)

for i = 2:N
    % 1. Update the heading direction with some noise
    hd(i) = hd(i-1) + randn() * turn_std;
    
    % 2. Calculate the proposed next position
    vx = speed * cos(hd(i));
    vy = speed * sin(hd(i));
    
    next_x = x(i-1) + vx * dt;
    next_y = y(i-1) + vy * dt;
    
    % 3. Wall Collision Logic (Bounce off the walls)
    % If it hits the Left or Right wall
    if next_x <= 0 || next_x >= boxSize
        hd(i) = pi - hd(i); % Reflect angle horizontally
        next_x = x(i-1) + (speed * cos(hd(i)) * dt); 
    end
    
    % If it hits the Top or Bottom wall
    if next_y <= 0 || next_y >= boxSize
        hd(i) = -hd(i);     % Reflect angle vertically
        next_y = y(i-1) + (speed * sin(hd(i)) * dt); 
    end
    
    % 4. Lock in the new position
    x(i) = next_x;
    y(i) = next_y;
end

% Ensure 'hd' is wrapped nicely between -pi and pi for the HD tuning curve later
hd = wrapToPi(hd);

%% 3. Calculate Instantaneous Firing Rates (Ground Truth)

% -- PLACE CELL RATE --
rate_PC = zeros(1, N);
for f = 1:size(pc_centers, 1)
    dist_sq = (x - pc_centers(f,1)).^2 + (y - pc_centers(f,2)).^2;
    rate_PC = rate_PC + exp(-dist_sq / (2 * pc_width^2));
end
rate_PC = (rate_PC / max(rate_PC)) * pc_max_rate;

% -- GRID CELL RATE (3-Cosine Interference Model) --
% Three waves offset by 60 degrees (pi/3)
angles = grid_orientation + [0, pi/3, 2*pi/3]; 
k = 4 * pi / (grid_spacing * sqrt(3)); % Wave number for spacing

rate_Grid = zeros(1, N);
for w = 1:3
    rate_Grid = rate_Grid + cos(k * (cos(angles(w)).*x + sin(angles(w)).*y));
end
% Exponentiate to create sharp distinct fields rather than blurry waves
rate_Grid = exp(rate_Grid);
rate_Grid = (rate_Grid / max(rate_Grid)) * grid_max_rate;

% -- HD CELL RATE (Von Mises Distribution) --
% Von Mises is the circular equivalent of a Gaussian
rate_HD = exp(hd_kappa * cos(hd - hd_pref));
rate_HD = (rate_HD / max(rate_HD)) * hd_max_rate;

%% 4. Inhomogeneous Poisson Spike Generator
% Probability of a spike in bin dt is roughly (rate * dt)
% Generate spikes where a random number is less than this probability
spikes_PC   = rand(1, N) < (rate_PC * dt);
spikes_Grid = rand(1, N) < (rate_Grid * dt);
spikes_HD   = rand(1, N) < (rate_HD * dt);

% Extract exact X, Y, and Time for the spikes
data_PC.t = t(spikes_PC); data_PC.x = x(spikes_PC); data_PC.y = y(spikes_PC);
data_Grid.t = t(spikes_Grid); data_Grid.x = x(spikes_Grid); data_Grid.y = y(spikes_Grid);
data_HD.t = t(spikes_HD); data_HD.x = x(spikes_HD); data_HD.y = y(spikes_HD);
data_HD.hd = hd(spikes_HD); % Store the angle the animal was facing

%% 5. Visualization
figure('Position', [100, 100, 1000, 300], 'Color', 'w');

% 1. Place Cell
ax1 = subplot(1, 3, 1);
plot(ax1, x, y, 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
hold(ax1, 'on');
scatter(ax1, data_PC.x, data_PC.y, 10, 'r', 'filled');
axis(ax1, 'square');
xlim(ax1, [0 boxSize]); 
ylim(ax1, [0 boxSize]);
title(ax1, sprintf('Place Cell (%d Spikes)', sum(spikes_PC)));

% 2. Grid Cell
ax2 = subplot(1, 3, 2);
plot(ax2, x, y, 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
hold(ax2, 'on');
scatter(ax2, data_Grid.x, data_Grid.y, 10, 'r', 'filled');
axis(ax2, 'square');
xlim(ax2, [0 boxSize]); 
ylim(ax2, [0 boxSize]);
title(ax2, sprintf('Grid Cell (%d Spikes)', sum(spikes_Grid)));

% 3. HD Cell (Polar Plot)
ax3 = subplot(1, 3, 3);
edges = linspace(-pi, pi, 30);
occupancy = histcounts(hd, edges) * dt;
spike_counts = histcounts(data_HD.hd, edges);
hd_tuning = spike_counts ./ (occupancy + 0.001); % Rate = Spikes / Time
centers = edges(1:end-1) + diff(edges)/2;

polarplot(centers, hd_tuning, 'LineWidth', 2, 'Color', 'b');
title(sprintf('HD Cell (%d Spikes)', sum(spikes_HD)));
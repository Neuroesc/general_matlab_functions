<img width="303" height="102" alt="neuroesc_long" src="https://github.com/user-attachments/assets/e185a933-e27d-4436-ab1f-52a4fe389e38" />

# Spatial cell simulations (Matlab)
A collection of MATLAB functions for simulating the firing properties of spatially modulated neurons, including place cells, grid cells, and head direction (HD) cells. This toolkit provides robust methods for generating synthetic data.

## Included Functions
### 1. Random Walk Generation
* **`gen_random_walk.m`**
    Simulates realistic animal exploratory trajectories inside arbitrary 2D geometries—from open fields to complex, multi-compartment mazes with walls and narrow doorways. Operating on a binary environment matrix (emap), it uses geodesic distance transforms and a weighted decision matrix to generate smooth, obstacle-aware movement paths. Key features:
1. Obstacle & doorway navigation: uses geodesic pathfinding (bwdistgeodesic) to select reachable macro-steps and interpolate collision-free paths around barriers and through narrow openings, this allows modelling complex geometries with small doorways.
2. Novelty-driven coverage: tracks recent spatial occupancy over a configurable memory window (novelty_min), biasing the agent away from recently visited areas to promote thorough exploration.
3. Tunable behavioral biases: combines seven weighted parameters—step distance, directional persistence (von Mises angular distribution), novelty, wall repulsion, central attraction, and custom repulsor/attractor maps—to shape movement dynamics.
4. Tracking-ready output: interpolates macro-steps to 50 Hz with optional positional jitter and smoothing (smoothn), returning an array of timestamps, smoothed/raw coordinates, head direction (degrees), and running speed (cm/s). 
<img width="1000" alt="image" src="https://github.com/user-attachments/assets/500bd7ba-0fa9-4315-9c5f-0ae731665222" />

<details>
  <summary>Click to view the code example</summary>

  ```matlab
% Open square
emap_square = false(122, 122);
emap_square(2:end-1, 2:end-1) = true;

% Square with dividing wall
emap_doorway = emap_square;
emap_doorway(:, 60:63) = false;      % 4 cm thick vertical dividing wall
emap_doorway(56:67, 60:63) = true;   % 12 cm doorway through the center

% Figure-8 maze
emap_fig8 = false(122, 122);
w = 14; % Corridor width in cm
mid = 61;
emap_fig8(2:(1+w), 2:end-1)             = true; % Top corridor
emap_fig8((end-w):end-1, 2:end-1)       = true; % Bottom corridor
emap_fig8(2:end-1, 2:(1+w))             = true; % Left corridor
emap_fig8(2:end-1, (end-w):end-1)       = true; % Right corridor
emap_fig8(2:end-1, (mid-w/2):(mid+w/2-1)) = true; % Central vertical stem

% Two connected rooms
emap_tworoom = false(82, 162);
emap_tworoom(2:81, 2:61)     = true; % Left room
emap_tworoom(2:81, 102:161)  = true; % Right room
emap_tworoom(35:48, 62:101)  = true; % 40 cm long, 14 cm wide corridor

% Get walk and plot
walk_mins = 10;
envs = {
    struct('name', 'Open Square',          'emap', emap_square,  'smoo', 10, 'center_pull', 100)
    struct('name', 'Divided Square (Door)', 'emap', emap_doorway, 'smoo', 4,  'center_pull', 0)
    struct('name', 'Figure-8 Corridors',   'emap', emap_fig8,    'smoo', 3,  'center_pull', 0)
    struct('name', 'Two Rooms + Corridor', 'emap', emap_tworoom, 'smoo', 4,  'center_pull', 0)
    };

f = figure('Name', 'Random Walk Geometries', 'Position', [100, 100, 1200, 900], 'Color', 'w');
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(envs)
    e = envs{i};
    disp(sprintf('\nSimulating %s', e.name));

    % Run random_walk (disable central attraction in multi-compartment mazes)
    [pos, ~] = gen_random_walk('emap', e.emap, ...
        'walk_length', walk_mins, ...
        'central_attraction', e.center_pull, ...
        'smoo', e.smoo, ...
        'novelty_min', 4);

    % Plot environment mask and trajectory
    % pos columns: [time, x_smooth, y_smooth, hd_deg, speed, x_raw, y_raw]
    nexttile;
    imagesc(e.emap);
    colormap(gca, [0.18 0.18 0.18; 0.98 0.98 0.98]); % Dark grey walls, white floor
    hold on;

    plot(pos(:, 2), pos(:, 3), 'Color', [0.05 0.45 0.85, 0.75], 'LineWidth', 1.1);
    plot(pos(1, 2), pos(1, 3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7);   % Start
    plot(pos(end, 2), pos(end, 3), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7); % End

    axis xy equal tight;
    set(gca, 'XTick', [], 'YTick', [], 'Box', 'on');
    title(sprintf('%s (%d min)', e.name, walk_mins), 'FontSize', 12, 'FontWeight', 'bold');
end
```
</details>

### 2. Synthetic Data Generation
* **`get_synth_neurons.m`**
  Generates a realistic correlated random walk trajectory for a virtual animal in a square arena. It also simulates synthetic spike trains for place cells, grid cells (using an interference model), and head direction Cells (using Von Mises tuning) via an inhomogeneous Poisson process. This is useful for unit-testing downstream analysis pipelines. It is a basic spatial simulation and does not model things like speed modulation, theta modulation or phase-precession.
<img width="1000" alt="image" src="https://github.com/user-attachments/assets/e77f169b-b887-4ba5-b055-3fc0af0ce24d" />

<details>
  <summary>Click to view the code example</summary>

  ```matlab
figure('Color','w')
tiledlayout(3,5,"TileSpacing",'tight')
N = 5;
[trajectory,cells] = get_synth_neurons(N,0,0,"hpc_widths",[2 40]);
for ii = 1:N
    [rmap,dmap] = rate_mapper([trajectory.x(:) trajectory.y(:)].*10,[cells.PC(ii).x(:) cells.PC(ii).y(:)].*10);
    m = get_spatial_info(dmap,rmap);
    nexttile
    imagesc(rmap)
    daspect([1 1 1])
    title(sprintf('Place Cell %d\nSpatial info: %.2fb/s',ii,m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
end

N = 5;
[trajectory,cells] = get_synth_neurons(0,N,0);
for ii = 1:N
    [rmap,dmap] = rate_mapper([trajectory.x(:) trajectory.y(:)].*10,[cells.Grid(ii).x(:) cells.Grid(ii).y(:)].*10);
    m = get_spatial_info(dmap,rmap);
    nexttile
    imagesc(rmap)
    daspect([1 1 1])
    title(sprintf('Grid Cell %d\nSpatial info: %.2fb/s',ii,m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
end

N = 5;
[trajectory,cells] = get_synth_neurons(0,0,N);
for ii = 1:N
    [rmap,dmap] = rate_mapper([trajectory.x(:) trajectory.y(:)].*10,[cells.HD(ii).x(:) cells.HD(ii).y(:)].*10);
    m = get_spatial_info(dmap,rmap);
    nexttile
    imagesc(rmap)
    daspect([1 1 1])
    title(sprintf('HD Cell %d\nSpatial info: %.2fb/s',ii,m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
end
```
</details>

## Requirements


## Basic Usage
Most functions are designed with an `arguments` block and support Name-Value pairs for flexible processing. Detailed usage instructions, expected inputs, and outputs can be found in the documentation block at the top of each file. 


function [trajectory,cells] = get_synth_neurons(nhpc,ngrid,nhd,opts)
% get_synth_neurons Simulates animal trajectory and spatial cells
%
% Generates a correlated random walk trajectory for a virtual animal in a 
% square arena, along with synthetic spike data for place cells, grid cells 
% (using an interference model), and head direction Cells (using Von Mises tuning). 
% Spike times are sampled via an inhomogeneous Poisson process based on the simulated trajectory.
%
% USAGE
%
% [trajectory, cells] = generate_synthetic_data(nhpc, ngrid, nhd) process with default settings
%
% [trajectory, cells] = generate_synthetic_data(__, Name, Value) process with Name-Value pairs 
%
% INPUT
%
% 'nhpc'      - Scalar, positive integer that specifies the number of Place Cells to generate.
%                Default value is 1.
%
% 'ngrid'    - Scalar, positive integer that specifies the number of Grid Cells to generate.
%                Default value is 1.
%
% 'nhd'      - Scalar, positive integer that specifies the number of Head Direction Cells to generate.
%                Default value is 1.
%
% 'duration'   - (Name-Value) Scalar, positive double that specifies the total session time.
%                Units are in seconds.
%                Default value is 600.
%
% 'boxsize'    - (Name-Value) Scalar, positive double that specifies the width/height of the square arena.
%                Units are in centimeters.
%                Default value is 100.
%
% 'fs'         - (Name-Value) Scalar, positive double that specifies the behavioral sampling rate.
%                Units are in Hertz (Hz).
%                Default value is 50.
%
% 'speed'      - (Name-Value) Scalar, positive double that specifies the animal's constant forward speed.
%                Units are in cm/s.
%                Default value is 20.
%
% 'plot'       - (Name-Value) Scalar logical that specifies whether to plot the first cell of each type.
%                Default value is false.
%
% 'hpc_fields' - 1x2, min and max number of place fields allowed
%                Default value is [1 2]
%
% 'hpc_widths' - 1x2, min and max place field width allowed, in cm
%                Default value is [8 18]
%
% 'grid_spacing' - 1x2, min and max grid field spacing allowed, in cm
%                Default value is [30 80]
%
% 'grid_max_angle' - 1x2, max grid orientation allowed, from x-axis, in radians
%                Default value is pi/3 (60 degrees)
%
% 'hd_widths' - 1x2, min and max head direction tuning width allowed
%                Default value is [2 6]
%
% OUTPUT
%
% 'trajectory' - Structure containing the fields 't', 'x', 'y', and 'hd'.
%                Units are in seconds, centimeters, and radians respectively.
%
% 'cells'      - Structure containing sub-structures 'PC', 'Grid', and 'HD'. Each contains 
%                the exact spike times ('t'), locations ('x', 'y'), and ground truth 'params' 
%                used to generate the synthetic neuron.
%
% NOTES
% 1. The trajectory is generated using a correlated random walk with perfectly 
%    elastic wall collisions to ensure uniform coverage of the environment.
%
% 2. Ground-truth parameters (like grid spacing or place field width) are drawn 
%    from random uniform distributions for each individual cell.
%
% EXAMPLE
% 
% % Generate 5 Place cells, 3 Grid cells, and 2 HD cells for a 20-minute session and plot them
% [traj, population] = generate_synthetic_data(5, 3, 2, 'duration', 1200, 'plot', true);
% 
% SEE ALSO 

% HISTORY
%
% version 1.0.0, Release 08/05/26 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        nhpc (1,1) double {mustBeNonnegative, mustBeInteger} = 1
        ngrid (1,1) double {mustBeNonnegative, mustBeInteger} = 1
        nhd (1,1) double {mustBeNonnegative, mustBeInteger} = 1
        
        % environment and simulation options
        opts.duration (1,1) double {mustBePositive} = 600 % sim duration in seconds
        opts.boxsize (1,1) double {mustBePositive} = 100  % environment size in cm
        opts.fs (1,1) double {mustBePositive} = 50        % sample rate of positions in Hz
        opts.speed (1,1) double {mustBePositive} = 20     % speed of the animal in cm/s
        opts.plot (1,1) logical = false                   % Plot first cell of each type?

        % cell limits and parameters
        opts.hpc_fields (1,2) double {mustBePositive} = [1 2] % min, max place fields
        opts.hpc_widths (1,2) double {mustBePositive} = [8 18] % min, max place field width
        opts.grid_spacing (1,2) double {mustBePositive} = [30 80] % min, max grid spacing
        opts.grid_max_angle (1,1) double {mustBePositive} = pi/3 % max grid orientation
        opts.hd_widths (1,2) double {mustBePositive} = [2 6] % min, max HD tuning width
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Generate Animal Trajectory (Correlated Random Walk)
    % Setup Time
    dt = 1 / opts.fs;
    t = 0:dt:opts.duration;
    N = length(t);
    boxSize = opts.boxsize;
    
    % preallocate
    x = zeros(1, N);
    y = zeros(1, N);
    hd = zeros(1, N);
    
    x(1) = boxSize / 2; y(1) = boxSize / 2; hd(1) = rand() * 2 * pi;
    turn_std = 0.2; 
    
    for i = 2:N
        hd(i) = hd(i-1) + randn() * turn_std;
        vx = opts.speed * cos(hd(i));
        vy = opts.speed * sin(hd(i));
        
        next_x = x(i-1) + vx * dt;
        next_y = y(i-1) + vy * dt;
        
        % Wall Collisions
        if next_x <= 0 || next_x >= boxSize
            hd(i) = pi - hd(i); 
            next_x = x(i-1) + (opts.speed * cos(hd(i)) * dt); 
        end
        if next_y <= 0 || next_y >= boxSize
            hd(i) = -hd(i);     
            next_y = y(i-1) + (opts.speed * sin(hd(i)) * dt); 
        end
        
        x(i) = next_x; y(i) = next_y;
    end
    hd = wrapToPi(hd);
    
    % Store trajectory output
    trajectory.t = t; 
    trajectory.x = x; 
    trajectory.y = y; 
    trajectory.hd = hd;
    
%%%%%%%%%%%%%%%% Generate place cells
    % Using struct arrays to hold N cells of each type
    cells.PC = struct('t', {}, 'x', {}, 'y', {}, 'params', {});
    cells.Grid = struct('t', {}, 'x', {}, 'y', {}, 'params', {});
    cells.HD = struct('t', {}, 'x', {}, 'y', {}, 'hd', {}, 'params', {});
    
    % place cells
    for i = 1:nhpc
        % random parameters
        num_fields = randi(opts.hpc_fields); % min or max fields
        centers = rand(num_fields, 2) * boxSize;
        width = opts.hpc_widths(1) + rand() * range(opts.hpc_widths); 
        max_rate = 10 + rand() * 20;      % 10 to 30 Hz
        
        % calculate rate
        rate = zeros(1, N);
        for f = 1:num_fields
            dist_sq = (x - centers(f,1)).^2 + (y - centers(f,2)).^2;
            rate = rate + exp(-dist_sq / (2 * width^2));
        end
        rate = (rate / max(rate)) * max_rate;
        
        % poisson spikes
        spikes = rand(1, N) < (rate * dt);
        
        % store
        cells.PC(i).t = t(spikes); 
        cells.PC(i).x = x(spikes); 
        cells.PC(i).y = y(spikes);
        cells.PC(i).params = struct('centers', centers, 'width', width, 'max_rate', max_rate, 'n_fields',num_fields);
    end
    
%%%%%%%%%%%%%%%% Generate grid cells
    for i = 1:ngrid
        % random parameters
        spacing = opts.grid_spacing(1) + rand() * range(opts.grid_spacing);
        orientation = rand() * opts.grid_max_angle;
        x0 = rand() * spacing;            % Random phase shift X
        y0 = rand() * spacing;            % Random phase shift Y
        max_rate = 10 + rand() * 15;      % 10 to 25 Hz
        
        % calculate rate
        angles = orientation + [0, pi/3, 2*pi/3]; 
        k = 4 * pi / (spacing * sqrt(3));
        rate = zeros(1, N);
        for w = 1:3
            rate = rate + cos(k * (cos(angles(w)).*(x-x0) + sin(angles(w)).*(y-y0)));
        end
        rate = exp(rate); % Sharpen fields
        rate = (rate / max(rate)) * max_rate;
        
        % poisson spikes
        spikes = rand(1, N) < (rate * dt);
        
        % store
        cells.Grid(i).t = t(spikes); 
        cells.Grid(i).x = x(spikes); 
        cells.Grid(i).y = y(spikes);
        cells.Grid(i).params = struct('spacing', spacing, 'orientation', orientation, 'max_rate', max_rate);
    end
    
%%%%%%%%%%%%%%%% Generate HD cells
    for i = 1:nhd
        % random parameters
        pref_dir = -pi + rand() * 2*pi;   % Anywhere from -pi to pi
        kappa = opts.hd_widths(1) + rand() * range(opts.hd_widths);
        max_rate = 15 + rand() * 25;      % 15 to 40 Hz
        
        % calculate rate (Von Mises)
        rate = exp(kappa * cos(hd - pref_dir));
        rate = (rate / max(rate)) * max_rate;
        
        % poisson spikes
        spikes = rand(1, N) < (rate * dt);
        
        % store
        cells.HD(i).t = t(spikes); 
        cells.HD(i).x = x(spikes); 
        cells.HD(i).y = y(spikes); 
        cells.HD(i).hd = hd(spikes);
        cells.HD(i).params = struct('pref_dir', pref_dir, 'kappa', kappa, 'max_rate', max_rate);
    end
    
%%%%%%%%%%%%%%%% Optional plotting
    if opts.plot
        figure('Position', [100, 100, 1000, 300], 'Color', 'w');
        
        if nhpc > 0
            ax1 = subplot(1, 3, 1);
            plot(ax1, x, y, 'Color', [0.8 0.8 0.8]); hold(ax1, 'on');
            scatter(ax1, cells.PC(1).x, cells.PC(1).y, 10, 'r', 'filled');
            axis(ax1, 'square'); 
            title(ax1, 'Place Cell 1');
        end
        
        if ngrid > 0
            ax2 = subplot(1, 3, 2);
            plot(ax2, x, y, 'Color', [0.8 0.8 0.8]); hold(ax2, 'on');
            scatter(ax2, cells.Grid(1).x, cells.Grid(1).y, 10, 'r', 'filled');
            axis(ax2, 'square'); 
            title(ax2, sprintf('Grid Cell 1 (Spacing: %.0f)', cells.Grid(1).params.spacing));
        end
        
        if nhd > 0
            subplot(1, 3, 3);
            edges = linspace(-pi, pi, 30);
            occupancy = histcounts(hd, edges) * dt;
            spike_counts = histcounts(cells.HD(1).hd, edges);
            tuning = spike_counts ./ (occupancy + 0.001);
            centers = edges(1:end-1) + diff(edges)/2;
            polarplot(centers, tuning, 'LineWidth', 2, 'Color', 'b');
            title('HD Cell 1');
        end
    end

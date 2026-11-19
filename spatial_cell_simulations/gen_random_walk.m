function trajectory = gen_random_walk(opts)
% gen_random_walk Generates a random walk within an arbitrary boundary matrix
%
% Models the exploratory trajectory of an animal navigating through simple or
% complex 2D environments. Uses geodesic distance transforms and a weighted
% decision matrix to route the agent around obstacles, through doorways, and
% away from recently visited areas (novelty preference).
%
% USAGE
%
% trajectory = gen_random_walk()
% trajectory = gen_random_walk(Name, Value)
%
% INPUT
%
% 'opts' - (Name-Value) Optional configuration parameters:
%   .emap               - (Logical matrix) Binary environment mask where true = walkable
%                         space and false = walls/obstacles.
%                         Default: 128x128 square (126x126 interior padded with 1-pixel wall).
%   .repulsor_map       - (Logical matrix) Mask of repulsor locations (same size as emap).
%                         Default: empty.
%   .attractor_map      - (Logical matrix) Mask of attractor locations (same size as emap).
%                         Default: empty.
%   .walk_length        - (Scalar double) Total duration of the walk in minutes.
%                         Default: 20.
%   .steps_per_min      - (Scalar double) Number of macro-decisions per minute.
%                         Default: 64.
%   .start_position     - (1x2 double) Initial [row, col] starting coordinates. If NaN or
%                         located on a wall, defaults to the nearest walkable pixel to center.
%                         Default: [NaN, NaN].
%   .start_angle        - (Scalar double) Initial heading angle in radians. If NaN, chosen
%                         uniformly at random from [-pi, pi].
%                         Default: NaN.
%   .step_dist          - (1x2 double) [mean, SD] of the Gaussian preferred macro-step
%                         distance (in pixels/cm).
%                         Default: [16, 8].
%   .step_fudge         - (1x2 double) [min, max] uniform random jitter added to the
%                         preferred step distance on each macro-step (in pixels/cm).
%                         Default: [-8, 8].
%   .angle_dist         - (1x2 double) [mean, kappa] of the von Mises distribution for
%                         turning angles relative to the previous heading (radians).
%                         Default: [0, 1].
%   .angle_fudge        - (1x2 double) [min, max] uniform random jitter added to the
%                         preferred turning angle on each macro-step (radians).
%                         Default: [-1.2, 1.2].
%   .smoo               - (Scalar double) Smoothing parameter passed to smoothn (robust).
%                         Set to 0 to disable smoothing.
%                         Default: 10.
%   .jit                - (Scalar double) Amplitude of uniform positional jitter added
%                         before smoothing to mimic tracking noise (0 to disable).
%                         Default: 0.
%   .central_attraction - (Scalar double) SD of Gaussian attraction toward the arena
%                         midpoint (smaller = stronger; 0 to disable).
%                         Default: 100.
%   .wall_repulsion     - (Scalar double) SD of Gaussian repulsion away from walls
%                         (smaller = stronger; 0 to disable).
%                         Default: 100.
%   .repulsion          - (Scalar double) SD of Gaussian repulsion from repulsor_map
%                         (smaller = stronger; 0 to disable).
%                         Default: 0.
%   .attraction         - (Scalar double) SD of Gaussian attraction to attractor_map
%                         (smaller = stronger; 0 to disable).
%                         Default: 0.
%   .novelty_min        - (Scalar double) Rolling memory window (minutes) used to track
%                         recent occupancy and bias exploration toward unvisited regions.
%                         Default: 8.
%   .decision_weights   - (1x7 double) Weights applied to each component of the decision
%                         matrix: [step_dist, angle, novelty, wall_repulsion,
%                         central_attraction, repulsors, attractors].
%                         Default: [1, 1, 1, 1, 1, 1, 1].
%
% OUTPUT
%
% 'trajectory' - (Struct) Simulated trajectory and session metadata with fields:
%   .t      - (Nx1 double) Timestamps in seconds.
%   .x      - (Nx1 double) Smoothed X coordinates (emap column space, pixels/cm).
%   .y      - (Nx1 double) Smoothed Y coordinates (emap row space, pixels/cm).
%   .v      - (Nx1 double) Instantaneous running speed (pixels/s or cm/s).
%   .raw_x  - (Nx1 double) Unsmoothed X coordinates prior to jitter/smoothing.
%   .raw_y  - (Nx1 double) Unsmoothed Y coordinates prior to jitter/smoothing.
%   .hd     - (Nx1 double) Movement heading direction in degrees [-180, 180]
%             (NaN at first and last samples).
%   .config - (Struct) Finalized configuration options used, including total .steps.
%   .srate  - (Scalar double) Trajectory sampling rate in Hz (50).
%
% EXAMPLE
%
% emap = padarray(ones(126, 'logical'), [1 1], false, 'both');
% traj = gen_random_walk('emap', emap);
%
% figure;
% tiledlayout('flow');
% nexttile;
% plot(traj.x, traj.y, 'k');
% daspect([1 1 1]);
% title('Square');
%
% nexttile;
% emap_door = padarray(ones(126, 'logical'), [1 1], false, 'both');
% emap_door(60, 1:60) = false;
% emap_door(60, 65:end) = false;
% traj_door = gen_random_walk('emap', emap_door);
% plot(traj_door.x, traj_door.y, 'k');
% daspect([1 1 1]);
% title('Narrow doorway');
%
% NOTES
% 1. Macro-steps are interpolated along shortest geodesic paths at 50 Hz.
% 2. X and Y coordinates map directly to columns and rows of emap (1 pixel = 1 cm).
% 3. Requires Image Processing Toolbox (bwdist, bwdistgeodesic, bwmorph, imgaussfilt),
%    Statistics and Machine Learning Toolbox (normpdf), and smoothn (when smoo > 0).
%
% SEE ALSO get_synth_neurons

% HISTORY
%
% version 2.2.0, Release 19/11/26 Renamed to gen_random_walk; switched output to trajectory struct
% version 2.1.0, Release 14/11/26 arguments block, optimized occupancy, meshgrid caching
% version 2.0.0, Release 14/11/26 Initial release, created from random_walk
% version 1.0.0, Release 08/10/24 random_walk created
%
% AUTHOR
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Memory and Place Systems (MAPS) lab
% roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        opts.emap (:,:) logical = padarray(ones(126,'logical'),[1 1],false,'both'); 
        opts.repulsor_map (:,:) logical = logical.empty;
        opts.attractor_map (:,:) logical = logical.empty;
        opts.walk_length (1,1) double = 20;  
        opts.steps_per_min (1,1) double = 64;  
        opts.start_position (1,2) double = [NaN NaN];
        opts.start_angle (1,1) double = NaN;    
        opts.step_dist (1,2) double = [16 8]; 
        opts.step_fudge (1,2) double = [-8 8]; 
        opts.angle_dist (1,2) double = [0 1]; 
        opts.angle_fudge (1,2) double = [-1.2 1.2]; 
        opts.smoo (1,1) double = 10; 
        opts.jit (1,1) double = 0; 
        opts.central_attraction (1,1) double = 100; 
        opts.wall_repulsion (1,1) double = 100; 
        opts.repulsion (1,1) double = 0; 
        opts.attraction (1,1) double = 0; 
        opts.novelty_min (1,1) double = 8; 
        opts.decision_weights (1,7) double = [1, 1, 1, 1, 1, 1, 1]; 
    end
    
    config = opts;
    config.steps = config.walk_length * config.steps_per_min;
    
    % Initialize maps
    wmap = config.emap;
    if isempty(config.repulsor_map); config.repulsor_map = zeros(size(wmap)); end
    if isempty(config.attractor_map); config.attractor_map = zeros(size(wmap)); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    fprintf('\tgenerating random walk (%.1f mins, %d steps)...\n', config.walk_length, config.steps);

    % Ready start point
    if all(isnan(config.start_position)) || ~wmap(config.start_position(1), config.start_position(2))
        lpoint = round(size(wmap) ./ 2);
        if ~wmap(lpoint(1), lpoint(2))
            cpoint = false(size(wmap));
            cpoint(lpoint(1), lpoint(2)) = true;
            d = bwdist(cpoint, "euclidean");
            d(~wmap) = NaN;
            [~, ind] = min(d, [], "all", "omitmissing", "linear");
            [r, c] = ind2sub(size(d), ind);
            lpoint = [r, c];
        end
    else
        lpoint = config.start_position;
    end

    if isnan(config.start_angle)    
        langle = -pi + (2*pi*rand()); 
    else
        langle = config.start_angle;
    end

    % Preallocate data structures
    steps = config.steps_per_min * config.walk_length;
    sample_rate = 50;
    samples_per_step = ceil((60 / config.steps_per_min) * sample_rate);
    pos = NaN(samples_per_step, steps, 3);
    
    step_lengths = config.step_dist(1) + (config.step_fudge(1) + (config.step_fudge(2)-config.step_fudge(1)).*rand(steps+1,1));
    step_angles  = config.angle_dist(1) + (config.angle_fudge(1) + (config.angle_fudge(2)-config.angle_fudge(1)).*rand(steps+1,1));

    % Precompute coordinate meshgrid for angle calculations (Massive Speedup)
    [cc, rr] = meshgrid(1:size(wmap,2), 1:size(wmap,1));

    % Precompute static environment repulsions/attractions
    if config.wall_repulsion
        ds = bwdist(~wmap, 'quasi-euclidean');
        wall_repulsion_map = 1 - rescale(normpdf(ds, 0, config.wall_repulsion), 'InputMin', 0); 
    else
        wall_repulsion_map = zeros(size(wmap));
    end

    if config.central_attraction
        wmap2 = false(size(wmap));
        mpoint = round(size(wmap)./2);
        wmap2(mpoint(1), mpoint(2)) = true;
        ds = bwdistgeodesic(wmap, wmap2, 'quasi-euclidean'); 
        central_attraction_map = rescale(normpdf(ds, 0, config.central_attraction), 'InputMin', 0); 
    else
        central_attraction_map = zeros(size(wmap));
    end

    repulsion_map = zeros(size(wmap));
    if config.repulsion && any(config.repulsor_map, "all")
        ds = bwdistgeodesic(wmap, config.repulsor_map, 'quasi-euclidean');
        repulsion_map = 1 - rescale(normpdf(ds, 0, config.repulsion), 'InputMin', 0); 
    end

    attraction_map = zeros(size(wmap));
    if config.attraction && any(config.attractor_map, "all")
        ds = bwdistgeodesic(wmap, config.attractor_map, 'quasi-euclidean'); 
        attraction_map = 1 - rescale(normpdf(ds, 0, config.attraction), 'InputMin', 0); 
    end

    % Memory optimization for occupancy window
    occ_counts = zeros(size(wmap));
    occ_history_steps = ceil(config.novelty_min * config.steps_per_min); 
    step_history_r = NaN(steps, 1);
    step_history_c = NaN(steps, 1);

    cpoint1 = false(size(wmap));  
    cpoint2 = false(size(wmap));            

    loopout = get_loop_progress(steps);    
    
    for ss = 1:steps
        cpoint1(:) = false; 
        cpoint1(lpoint(1), lpoint(2)) = true; 

        % Distance logic
        ds1 = bwdistgeodesic(wmap, cpoint1, 'quasi-euclidean');
        ds_map = rescale(normpdf(ds1, step_lengths(ss), config.step_dist(2)), 'InputMin', 0); 

        % Angle logic (using native wrapToPi instead of angdiff)
        as1 = atan2(rr - lpoint(1), cc - lpoint(2)); 
        as2 = wrapToPi(as1 - langle); 
        as_map = rescale(circ_vmpdf(as2, step_angles(ss), config.angle_dist(2)), 'InputMin', 0); 

        % Optimized macro-step occupancy logic
        occ = imgaussfilt(occ_counts, 5);       
        occ_map = 1 - rescale(occ, 'InputMin', 0); 

        % Combine decision matrix
        dm = (config.decision_weights(1).*ds_map) + ...
             (config.decision_weights(2).*as_map) + ...
             (config.decision_weights(3).*occ_map) + ...
             (config.decision_weights(4).*wall_repulsion_map) + ...
             (config.decision_weights(5).*central_attraction_map) + ...
             (config.decision_weights(6).*repulsion_map) + ...
             (config.decision_weights(7).*attraction_map);
         
        dm(isinf(ds1) | ~wmap) = NaN;
        dm(lpoint(1), lpoint(2)) = NaN;
        
        [~, ind] = max(dm, [], "all", "omitmissing", "linear");
        [r, c] = ind2sub(size(dm), ind);
        lpoint_next = [r, c];
        langle_next = as1(ind);

        % Update occupancy history efficiently
        step_history_r(ss) = r;
        step_history_c(ss) = c;
        occ_counts(r, c) = occ_counts(r, c) + 1;
        
        old_idx = ss - occ_history_steps;
        if old_idx > 0 && ~isnan(step_history_r(old_idx))
            old_r = step_history_r(old_idx);
            old_c = step_history_c(old_idx);
            occ_counts(old_r, old_c) = occ_counts(old_r, old_c) - 1;
        end

        % Path logic
        cpoint2(:) = false; 
        cpoint2(r, c) = true; 
        ds2 = bwdistgeodesic(wmap, cpoint2, 'quasi-euclidean');   
        ds3 = round(ds1 + ds2); 
        ds3(isnan(ds3)) = max(ds3, [], "all", "omitmissing");
        
        paths = imregionalmin(ds3); 
        paths_thinned = bwmorph(paths, 'thin', inf); 
        paths_thinned(lpoint_next(1), lpoint_next(2)) = true; 
        paths_thinned(lpoint(1), lpoint(2)) = true; 

        vals = [rr(paths_thinned(:)), cc(paths_thinned(:)), ds1(paths_thinned(:))];
        [~, sidx] = sort(vals(:,3));
        vals = vals(sidx, :);

        ipath = interp1(1:size(vals,1), vals(:,1:2), linspace(1, size(vals,1), samples_per_step+1), 'linear');
        ipath = ipath(2:end, :); 
        
        t = (1:size(ipath,1)) .* (1/sample_rate);
        if ss > 1
            t = t + max(pos(:, ss-1, 1), [], "all", "omitmissing");
        end
        pos(:, ss, :) = [t(:), ipath];

        lpoint = lpoint_next;
        langle = langle_next;
        loopout = get_loop_progress(loopout);    
    end

    % Smooth and construct final output array
    pos_unwrapped = reshape(pos(:), [], 3, 1);
    pox0 = pos_unwrapped(:,3); 
    poy0 = pos_unwrapped(:,2); 

    poxn = pox0 + (rand(size(pox0))-0.5) .* config.jit; 
    poyn = poy0 + (rand(size(poy0))-0.5) .* config.jit; 

    if config.smoo
        d = smoothn({poxn, poyn}, config.smoo, 'robust');
        pox = d{1}; poy = d{2};  
    else
        pox = poxn; poy = poyn;
    end
    pot = pos_unwrapped(:,1);

    % Heading and speed
    p_now = [pox(:), poy(:)];
    p_next = circshift(p_now, -1, 1); 
    p_last = circshift(p_now, 1, 1); 
    
    [poh, ds] = cart2pol(p_next(:,1) - p_last(:,1), p_next(:,2) - p_last(:,2)); 
    poh([1 end]) = NaN;
    poh = rad2deg(poh);  

    wsize = 3;
    ds([1 end]) = NaN;
    ds = ds ./ 2; 
    dist_sum = movsum(ds, wsize, 'omitnan', 'Endpoints', 'shrink'); 
    pov = dist_sum ./ (wsize * (1/sample_rate)); 

    % Save output
    trajectory = struct;
    trajectory.t = pot(:); 
    trajectory.x = pox(:); 
    trajectory.y = poy(:); 
    trajectory.v = pov(:); 
    trajectory.raw_x = pox0(:); 
    trajectory.raw_y = poy0(:);     
    trajectory.hd = poh(:); % in degrees
    trajectory.config = config;
    trajectory.srate = sample_rate;

end

% Von Mises Distribution for native angle distribution
function p = circ_vmpdf(alpha, thetahat, kappa)
    C = 1 ./ (2 * pi * besseli(0, kappa));
    p = C * exp(kappa * cos(alpha - thetahat));
end
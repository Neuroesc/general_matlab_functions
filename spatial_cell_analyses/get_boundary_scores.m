function [allo_stats, ego_stats, speedlift, maps] = get_boundary_scores(pos, sindx, epoly, sample_rate, opts)
% calc_boundary_ratemaps Computes allocentric and egocentric boundary tuning
%
% This function constructs boundary ratemaps by calculating the distance 
% and angle to environmental borders for every position frame and spike. 
% It computes both allocentric and egocentric tuning metrics including Mean 
% Resultant Length (MRL), preferred angle, and preferred distance. It utilizes 
% a speedlift structure to cache computationally expensive raycasting 
% operations per session to accelerate multi-cell processing.
%
% USAGE
%
% [allo_stats, ego_stats, speedlift] = get_boundary_scores(pos, sindx, epoly, 50)
% [allo_stats, ego_stats, speedlift] = get_boundary_scores(pos, sindx, epoly, 50, 'speedlift', speedlift_data, 'max_dist', 50)
%
% INPUT
%
% 'pos'         - (numeric) N x 3 matrix containing [X, Y, Head Direction (rad)].
% 'sindx'     - (numeric) Indices of 'pos' where spikes occurred.
% 'epoly'       - (numeric) K x 2 matrix of environment boundary vertices.
% 'sample_rate' - (numeric) Position sampling rate in Hz (e.g., 50).
%
% OPTIONAL NAME-VALUE INPUTS
%
% 'speedlift'    - (struct) Cached distance matrices. Pass empty to generate.
% 'max_dist'     - (numeric) Maximum boundary distance to compute in cm (default: 62.5).
% 'ang_step'     - (numeric) Angular bin size in degrees (default: 3).
% 'dist_step'    - (numeric) Distance bin size in cm (default: 2.5).
% 'smooth_sigma' - (numeric) Gaussian smoothing standard deviation in bins (default: 5).
%
% OUTPUT
%
% 'allo_stats'  - (struct) Allocentric boundary ratemap, MRL, preferred angle and distance.
% 'ego_stats'   - (struct) Egocentric boundary ratemap, MRL, preferred angle and distance.
% 'speedlift'   - (struct) Cached distance and angle matrices for subsequent cells.
%
% NOTES
% 1. The Gaussian smoothing explicitly uses circular padding for the angular 
%    dimension and replicate padding for the distance dimension before applying 
%    imgaussfilt to prevent erroneous boundary bleeding.
%
% EXAMPLE
% 
% % 1. Run the function
% [allo, ego, speedlift] = get_boundary_scores(pos, sindx, epoly, 50, 'max_dist', 62.5);
% 
% % 2. Plot the Egocentric Ratemap in polar coordinates
% theta = deg2rad(0 : 3 : 357);
% r = 1.25 : 2.5 : 62.5; % Bin centers
% [R, Theta] = meshgrid(r, theta);
% [X, Y] = pol2cart(Theta, R);
% 
% figure;
% surf(X, Y, ego.ratemap, 'EdgeColor', 'none');
% view(2); axis equal; axis off; colormap('turbo');
% title('Egocentric Boundary Map');
% 
% % 3. Add Orientation Labels
% max_r = max(r) * 1.15;
% text(0, max_r, 'Front (0°)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% text(0, -max_r, 'Behind (180°)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% text(-max_r, 0, 'Left (90°)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% text(max_r, 0, 'Right (270°)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% 
% SEE ALSO calc_glm_predictors run_rsctp_pipeline

% HISTORY
%
% version 1.1.0, Release 25/08/26 Added tunable parameters and imgaussfilt support
% version 1.0.0, Release 25/08/26 Initial release
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        pos (:,3) double
        sindx (:,1) double
        epoly (:,2) double
        sample_rate (1,1) double
        
        opts.speedlift struct = struct()
        opts.max_dist (1,1) double = 62.5
        opts.ang_step (1,1) double = 3
        opts.dist_step (1,1) double = 2.5
        opts.smooth_sigma (1,1) double = 5
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % 1. Define binning parameters
    angles_deg = 0 : opts.ang_step : (360 - opts.ang_step);
    angles_rad = deg2rad(angles_deg);
    n_ang = length(angles_rad);
    
    ang_edges = deg2rad(0 : opts.ang_step : 360);
    dist_edges = 0 : opts.dist_step : opts.max_dist;
    
    n_frames = size(pos, 1);
    
    % 2. Speedlift Generation (Compute rays if not provided)
    if isempty(fieldnames(opts.speedlift)) || isempty(opts.speedlift)
        allo_dists = zeros(n_frames, n_ang);
        ego_dists  = zeros(n_frames, n_ang);
        
        % Convert HD to angular bins to shift allocentric rays to egocentric
        hd_rad = wrapTo2Pi(pos(:, 3));
        hd_bins = round(rad2deg(hd_rad) / opts.ang_step) + 1;
        hd_bins(hd_bins == (n_ang + 1)) = 1; 
        
        for t = 1:n_frames
            pos_t = pos(t, 1:2);
            dists = cast_rays(pos_t, epoly, angles_rad);
            allo_dists(t, :) = dists;
            
            % Shift allocentric rays by animal's HD to get egocentric rays
            % In this egocentric frame, 0/360 is front, 90 is left, 180 is behind
            shift_val = -(hd_bins(t) - 1); 
            ego_dists(t, :) = circshift(dists, shift_val);
        end
        
        % Cache for future cells
        speedlift.allo_dists = allo_dists;
        speedlift.ego_dists = ego_dists;
        speedlift.angles_mat = repmat(angles_rad, n_frames, 1);
    else
        speedlift = opts.speedlift;
    end
    
    % 3. Extract data from speedlift
    allo_dists = speedlift.allo_dists;
    ego_dists  = speedlift.ego_dists;
    angles_mat = speedlift.angles_mat;
    
    % 4. Compute Occupancy Maps (seconds per bin)
    [occ_allo, ~, ~] = histcounts2(angles_mat(:), allo_dists(:), ang_edges, dist_edges);
    [occ_ego, ~, ~] = histcounts2(angles_mat(:), ego_dists(:), ang_edges, dist_edges);
    
    occ_allo = occ_allo / sample_rate;
    occ_ego  = occ_ego / sample_rate;
    
    % 5. Compute Spike Maps
    spt_angles_mat = repmat(angles_rad, length(sindx), 1);
    
    spt_allo_dists = allo_dists(sindx, :);
    [spk_allo, ~, ~] = histcounts2(spt_angles_mat(:), spt_allo_dists(:), ang_edges, dist_edges);
    
    spt_ego_dists = ego_dists(sindx, :);
    [spk_ego, ~, ~] = histcounts2(spt_angles_mat(:), spt_ego_dists(:), ang_edges, dist_edges);
    
    % 6. Compute and Smooth Ratemaps
    raw_allo_map = spk_allo ./ max(occ_allo, eps);
    raw_ego_map  = spk_ego ./ max(occ_ego, eps);
    
    % Remove unvisited bins
    raw_allo_map(occ_allo < 0.05) = NaN;
    raw_ego_map(occ_ego < 0.05) = NaN;
    
    % Smooth maps using custom imgaussfilt wrapper
    smoothed_allo_rmap = smooth_boundary_map(raw_allo_map, opts.smooth_sigma);
    smoothed_allo_dmap = smooth_boundary_map(occ_allo, opts.smooth_sigma);
    smoothed_ego_rmap  = smooth_boundary_map(raw_ego_map, opts.smooth_sigma);
    smoothed_ego_dmap  = smooth_boundary_map(occ_ego, opts.smooth_sigma);

    % 7. Calculate Tuning Metrics
    dist_centers = dist_edges(1:end-1) + (opts.dist_step / 2);
    allo_stats = compute_tuning_metrics(smoothed_allo_rmap, smoothed_allo_dmap, angles_rad, dist_centers);
    ego_stats  = compute_tuning_metrics(smoothed_ego_rmap, smoothed_ego_dmap, angles_rad, dist_centers);
    
    allo_stats.ratemap = smoothed_allo_rmap;
    allo_stats.dwellmap = smoothed_allo_dmap;
    ego_stats.ratemap = smoothed_ego_rmap;
    ego_stats.dwellmap = smoothed_ego_dmap;    
    
end

%%%%%%%%%%%%%%%% HELPER FUNCTIONS
function dists = cast_rays(pos, epoly, angles)
    % Casts rays from position to environmental polygon boundaries
    n_rays = length(angles);
    dists = inf(1, n_rays);
    for i = 1:size(epoly,1)-1
        p1 = epoly(i,:);
        p2 = epoly(i+1,:);
        v1 = pos - p1;
        v2 = p2 - p1;
        v3 = [sin(angles)', -cos(angles)']; 
        dot_val = v2(1)*v3(:,1) + v2(2)*v3(:,2);
        t1 = (v2(1)*v1(2) - v2(2)*v1(1)) ./ dot_val; 
        t2 = (v1(1)*v3(:,1) + v1(2)*v3(:,2)) ./ dot_val; 
        hit_idx = t1 > 0 & t2 >= 0 & t2 <= 1;
        dists(hit_idx) = min(dists(hit_idx), t1(hit_idx)');
    end
end

function smoothed = smooth_boundary_map(raw_map, sigma)
    % Fills NaNs, applies mixed-padding, and smooths with a 2D Gaussian
    
    % Temporarily fill NaNs with 0 for convolution
    filled_map = raw_map;
    filled_map(isnan(filled_map)) = 0;
    
    % Define kernel and padding size based on standard deviation
    pad_sz = ceil(3 * sigma);
    kernel_sz = 2 * pad_sz + 1;
    H = fspecial('gaussian', [kernel_sz kernel_sz], sigma);
    
    % Pad circularly on dimension 1 (Angles) and replicate on dimension 2 (Distance)
    padded = padarray(filled_map, [pad_sz 0], 'circular', 'both');
    padded = padarray(padded, [0 pad_sz], 'replicate', 'both');
    
    % Apply 2D convolution, returning only the valid (unpadded) region
    smoothed = conv2(padded, H, 'valid');
    
    % Re-apply NaNs to unvisited bins
    smoothed(isnan(raw_map)) = NaN;
end

function stats = compute_tuning_metrics(rmap, dmap, angles_rad, dist_centers)
    % Computes MRL, MRA, and preferred distance from the smoothed ratemap
    
    % spatial information content etc
    stats = get_spatial_info(dmap,rmap,'metrics',{'spatial_info','kld'});

    % Collapse across distance to get angular tuning vector
    F_theta = mean(rmap, 2, 'omitnan');
    F_theta(isnan(F_theta)) = 0;
    
    % Calculate Mean Resultant (MR) and MRL
    MR = sum(F_theta .* exp(1i * angles_rad')) / sum(F_theta);
    stats.MRL = abs(MR);
    
    % Calculate Mean Resultant Angle (Preferred Angle)
    stats.pref_angle_rad = wrapTo2Pi(angle(MR));
    stats.pref_angle_deg = rad2deg(stats.pref_angle_rad);
    
    % Find preferred distance at the closest angular bin to the MRA
    [~, mra_idx] = min(abs(angles_rad - stats.pref_angle_rad));
    
    % Extract distance tuning curve at that specific angle
    dist_tuning = rmap(mra_idx, :);
    [~, max_dist_idx] = max(dist_tuning);
    
    if ~isempty(max_dist_idx) && ~isnan(dist_tuning(max_dist_idx))
        stats.pref_dist = dist_centers(max_dist_idx);
    else
        stats.pref_dist = NaN;
    end
end
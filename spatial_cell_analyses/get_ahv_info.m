function [scores, slopes, ahv_results] = get_ahv_info(pot, poh, spindx, opts)
% calc_ahv_tuning Calculates angular head velocity maps and tuning scores
%
% This function calculates angular head velocity (AHV), constructs AHV 
% dwell time maps, spike maps, and firing rate maps, and computes AHV 
% tuning scores and slopes. AHV is calculated as the first derivative 
% over a 200-ms sliding window ending on the current time point. Positive 
% velocities denote right turns, and negative velocities denote left turns.
%
% USAGE
%
% [scores, slopes, res] = get_ahv_info(pot, poh, spindx, opts);
%
% INPUT
%
% 'pot'    - (double) Nx1 array of position time stamps in seconds.
%
% 'poh'    - (double) Nx1 array of head direction angles in degrees 
%                 (-180 to 180, CCW is positive).
%
% 'spindx' - (double) Sx1 spike index into position data, i.e. the nearest
%           neighbour in pot for each spike
%
% 'max_ahv' - (double) Maximum AHV in tuning curves (deg/s)
%               default = 80
%
% 'bin_size' - (double) Binsize for tuning curves (deg/s)
%               default = 6
%
% 'window_ms' - (double) Calculate first derivative over this time window (ms)
%               default = 200
%
% 'min_dwell' - (double) Bins with less than this dwell time will be considered
%               unvisited (seconds)
%               default = 0.05
%
% OUTPUT
%
% 'scores'  - (double) 1x2 array, left AHV score (Pearson's r) and right AHV
%               score
%
% 'slopes'  - (double) 1x2 array, left AHV slope (linear fit) and right AHV
%               slope
%
% 'ahv_results' - (struct) A structure containing the dwell map, spike map, 
%                 rate map, bin centers, and left/right Pearson's r 
%                 correlations and linear slopes.
%
% NOTES
% 1. Method follows Keshavarzi et al. (2021). AHV tuning curves use a 
%    default maximum angular velocity of 80 degrees/s binned at 6 degrees/s.
%
% 2. Assumes an approximately uniform sampling rate for pot (e.g., 50Hz).
%
% EXAMPLE
% 
% % Calculate AHV tuning with defaults
% [scores, slopes, res] = get_ahv_info(pot, poh, spindx);
% 
% SEE ALSO get_directional_info

% HISTORY
%
% version 1.0.0, Release 24/08/26 Initial release
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        pot (:,1) double
        poh (:,1) double
        spindx (:,1) double
        opts.max_ahv (1,1) double = 80
        opts.bin_size (1,1) double = 6
        opts.window_ms (1,1) double = 200
        opts.min_dwell (1,1) double = 0.05        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Calculate Angular Head Velocity (AHV)
    % Determine sampling interval (dt) and window size in samples
    dt = median(diff(pot));
    window_samples = round((opts.window_ms / 1000) / dt);
    
    % Unwrap head direction to prevent artifactual jumps at -180/180
    hd_unwrapped = unwrap(deg2rad(poh));
    
    % Calculate the derivative over a 200-ms window ending on that time point
    ahv_rad = zeros(length(hd_unwrapped), 1);
    ahv_rad(window_samples+1:end) = (hd_unwrapped(window_samples+1:end) - hd_unwrapped(1:end-window_samples)) / (opts.window_ms / 1000);
    
    ahv_deg = rad2deg(ahv_rad);
    
    % Convert sign: Original data is CCW positive (Left turn = positive). 
    % Keshavarzi et al. assigned positive velocities to right turns.
    ahv_deg = -ahv_deg; 
    
    % 2. Define Binning Edges
    % Compute edges safely encompassing max_ahv using the specified bin_size
    edge_limit = ceil(opts.max_ahv / opts.bin_size) * opts.bin_size;
    edges = -edge_limit : opts.bin_size : edge_limit;
    bin_centers = edges(1:end-1) + (opts.bin_size / 2);
    
    % 3. Construct Maps
    % Dwell time map
    dwell_map = histcounts(ahv_deg, edges) .* dt;
    
    % % Assign spikes to nearest position timestamps
    % spike_idx = interp1(pot, 1:length(pot), spt, 'nearest');
    % spike_idx(isnan(spike_idx)) = []; % Clean out-of-bounds spikes
    
    % Spike map and Firing rate map
    spike_map = histcounts(ahv_deg(spindx), edges);
    rate_map = spike_map ./ dwell_map;
    
    % Handle unvisited bins (prevent NaNs)
    rate_map(dwell_map < opts.min_dwell) = NaN; % threshold unvisited bins (< 50ms)
    
    % 4. Compute AHV Scores and Slopes
    % Split into left (negative) and right (positive) turns
    right_idx = bin_centers > 0;
    left_idx = bin_centers < 0;
    
    % Helper inline function for safe linear fit and correlation
    [right_r, right_slope] = fit_tuning_curve(bin_centers(right_idx), rate_map(right_idx));
    [left_r, left_slope]   = fit_tuning_curve(bin_centers(left_idx), rate_map(left_idx));
    
    % spatial information content etc
    m = get_spatial_info(dwell_map,rate_map,'metrics',{'spatial_info','kld'});

    % 5. Store Results
    ahv_results = m;
    ahv_results.edges = edges;
    ahv_results.bin_centers = bin_centers;
    ahv_results.dwell_map = dwell_map;
    ahv_results.spike_map = spike_map;
    ahv_results.rate_map = rate_map;
    ahv_results.right_score = right_r;
    ahv_results.left_score = left_r;
    ahv_results.right_slope = right_slope;
    ahv_results.left_slope = left_slope;
    scores = [left_r right_r];
    slopes = [left_slope right_slope];
end

%%%%%%%%%%%%%%%% HELPER FUNCTIONS
function [abs_r, slope] = fit_tuning_curve(x, y)
    % Removes NaNs and computes absolute Pearson's r and linear slope
    valid = ~isnan(x) & ~isnan(y) & ~isinf(y);
    x_val = x(valid);
    y_val = y(valid);
    
    if length(x_val) < 3
        abs_r = NaN;
        slope = NaN;
        return;
    end
    
    % Absolute Pearson's correlation (magnitude)
    r_mat = corrcoef(x_val, y_val);
    abs_r = abs(r_mat(1,2));
    
    % Linear slope via polyfit
    p = polyfit(x_val, y_val, 1);
    slope = p(1);
end
function [score, slope, intercept, speed_results] = get_speed_info(pot, pox, poy, spindx, opts)
% get_speed_info Calculates speed tuning metrics according to Kropff et al. (2015)
%
% This function calculates instantaneous running speed, applies the necessary 
% temporal smoothing, and constructs the required rate maps based on Kropff 
% et al. (2015). It computes the speed tuning curve, speed score (Pearson 
% correlation), and linear regression parameters (slope and intercept).
%
% USAGE
%
% [score, slope, intercept, res] = get_speed_info(pot, pox, poy, spindx)
% [score, slope, intercept, res] = get_speed_info(pot, pox, poy, spindx, opts)
%
% INPUT
%
% 'pot'           - (double) Nx1 array of position time stamps in seconds.
%
% 'pox'         - (double) Nx1 array of X position coordinates in cm.
%
% 'poy'         - (double) Nx1 array of Y position coordinates in cm.
%
% 'spindx' - (uint32/double) Sx1 array of spike indices into position data.
%
% 'bin_size' - (double) Bin size to use for speed tuning curve (cm/s)
%               default = 2
%
% 'max_speed' - (double) Maximum speed to calculate tuning curve to (cm/s)
%               default = 50
%
% 'smooth_sigma_speed' - (double) Smooth speed maps by this amount
%               Gaussian (sd)
%               default = 3
%
% 'smooth_width_fr' - (double) Smooth instantaneous firing rate by this amount
%               Gaussian width (s)
%               default = 0.250
%
% 'min_dwell_prc' - (double) Bins with less than this % of total dwell time will
%               be considered unvisited (%)
%               default = 0.5
%
% OUTPUT
%
% 'score'         - (double) Speed score (Pearson product-moment correlation).
% 'slope'         - (double) Slope of the speed tuning curve linear fit.
% 'intercept'     - (double) Intercept of the speed tuning curve linear fit.
% 'speed_results' - (struct) Structure containing tuning curves and maps.
%
% NOTES
% 1. Method follows Kropff et al. (2015). Tuning curves use equally spaced 
%    bins (2 cm/s) smoothed by a Gaussian filter (SD: 3 cm/s).
%
% 2. Bins with less than 0.5% data coverage are excluded as valid.
%
% 3. Static periods (< 2 cm/s) are filtered out prior to calculating the score 
%    to avoid behavioral state confounds in open field recordings.
%
% EXAMPLE
% 
% % Calculate speed tuning metrics with defaults
% [score, slope, intercept, res] = get_speed_info(pot, pox, poy, spindx);
% 
% SEE ALSO get_ahv_info

% HISTORY
%
% version 1.0.1, Release 24/08/26 Switched to spt_pot_index for direct 50Hz frame alignment
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
        pox (:,1) double
        poy (:,1) double
        spindx (:,1) double
        opts.bin_size (1,1) double = 2 % cm/s
        opts.max_speed (1,1) double = 50 % cm/s
        opts.smooth_sigma_speed (1,1) double = 3 % cm/s
        opts.smooth_width_fr (1,1) double = 0.250 % 250 ms
        opts.min_dwell_prc (1,1) double = 0.5 % 0.5% coverage
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % 1. Calculate Instantaneous Speed
    dt = median(diff(pot));
    dx = [0; diff(pox)];
    dy = [0; diff(poy)];
    inst_speed = sqrt(dx.^2 + dy.^2) ./ dt; 

    % 2. Calculate Instantaneous Firing Rate (50 Hz frames / 20 ms bins)
    % Filter out-of-bound indices if any
    valid_spt_idx = spindx(spindx >= 1 & spindx <= length(pot));
    spike_counts = accumarray(valid_spt_idx, 1, [length(pot), 1]);
    
    % Smooth with 250-ms-wide Gaussian filter (Kropff et al. 2015)
    window_len = round(opts.smooth_width_fr / dt);
    gauss_win = gausswin(window_len);
    gauss_win = gauss_win / sum(gauss_win);
    inst_fr = conv(spike_counts ./ dt, gauss_win, 'same');

    % 3. Speed Score (Pearson correlation)
    % Filter out static periods (< 2 cm/s) for open field recordings
    % From Kropff et al. (2015):
    % The speed score for each cell was defined as the Pearson product-moment correlation 
    % between the cell's instantaneous firing rate and the rat's instantaneous running 
    % speed, on a scale from −1 to 1.
    % ... and ...
    % In order to distinguish speed-correlated effects from changes in behavioural state 
    % (foraging versus sitting still), we dismissed all data produced at a running speed 
    % lower than 2 cm s−1 in the calculation of the observed and shuffled speed scores.    
    run_idx = inst_speed >= 2 & inst_speed < opts.max_speed;
    if sum(run_idx) > 2
        score = corr(inst_speed(run_idx), inst_fr(run_idx),"Type","Pearson","Rows","pairwise");
    else
        score = NaN;
    end

    % figure
    % plot(inst_speed(run_idx),inst_fr(run_idx),'ko')
    % keyboard

    % 4. Construct Speed Rate Map
    % From Kropff et al. (2015):
    % Rate maps that showed firing rate as a function of location, head direction or speed 
    % were constructed with similar procedures. Histograms of spike count on one hand and 
    % time spent in the location on the other were built for each cell, using equally spaced 
    % bins (bin size: 2.5 cm for spatial maps, 6° for head-direction maps and 2 cm s−1 for 
    % speed maps). Each bin of the rate map was obtained as the ratio between the spike count 
    % and the time spent, smoothed by a Gaussian filter (standard deviation: 4 cm for spatial 
    % maps, 6° for head-direction maps and 3 cm s−1 for speed maps). In speed maps, where 
    % coverage is very inhomogeneous, only bins accounting for at least 0.5% of the 
    % data were included as valid. In composite rate maps for speed versus head direction, 
    % this threshold was divided by the number of head-direction bins. Instantaneous firing 
    % rate was obtained by dividing the whole session into 20-ms bins, coinciding with the 
    % frames of the tracking camera. A temporal histogram of spiking was then obtained, 
    % smoothed with a 250-ms-wide Gaussian filter. Spatial and head-directional information 
    % measures32 were based on these maps.
    edges = 0:opts.bin_size:opts.max_speed;
    bin_centers = edges(1:end-1) + (opts.bin_size / 2);
    
    [dwell_counts, ~, bin_idx] = histcounts(inst_speed, edges);
    dwell_map = dwell_counts * dt;
    
    % Accumulate frame spike counts into speed bins
    spike_map = accumarray(bin_idx(bin_idx > 0), spike_counts(bin_idx > 0), [length(bin_centers) 1])';
    raw_rate_map = spike_map ./ dwell_map;

    % Gaussian smoothing for speed maps (SD: 3 cm/s)
    sigma_bins = opts.smooth_sigma_speed / opts.bin_size;
    rate_map = imgaussfilt(raw_rate_map, sigma_bins, 'Padding', 'replicate');

    % Invalidate bins with < 0.5% data coverage
    valid_bins = (dwell_map ./ sum(dwell_map)) >= (opts.min_dwell_prc/100);
    rate_map(~valid_bins) = NaN;

    % 5. Linear Regression (Slope and Intercept)
    valid_x = bin_centers(valid_bins);
    valid_y = rate_map(valid_bins);
    if length(valid_x) >= 2
        p = polyfit(valid_x, valid_y, 1);
        slope = p(1);
        intercept = p(2);
    else
        slope = NaN;
        intercept = NaN;
    end

    % spatial information content etc
    m = get_spatial_info(dwell_map,rate_map,'metrics',{'spatial_info','kld'});

    % 6. Package Results
    speed_results = m;
    speed_results.speed_score = score;    
    speed_results.slope = slope;    
    speed_results.y_intercept = intercept;    
    speed_results.inst_speed = inst_speed;
    speed_results.inst_fr = inst_fr;
    speed_results.bin_centers = bin_centers;
    speed_results.rate_map = rate_map;
    speed_results.dwell_map = dwell_map;
end
function [res,dat] = get_spatial_bootshuff(pox,poy,pot,rmset,opts)
% get_spatial_bootshuff Computes bootstrapped and shuffled spatial metrics
%
% Calculates robust spatial, grid, and directional metrics for a single cell. 
% It uses bootstrapping (spike resampling with replacement) to determine stable 
% observed values, and temporal shuffling (circular time shifts) to generate 
% a null distribution. It outputs both parametric and non-parametric significance.
% This is based on the Method described by Savelli et al. (2017), 
% https://elifesciences.org/articles/21354, with some improvements.
%
% USAGE
%
% [res, dat] = get_spatial_bootshuff(pox, poy, pot, rmset) process with default settings
%
% [res, dat] = get_spatial_bootshuff(pox, poy, pot, rmset, 'sindx', spike_indices) process using pre-computed indices
%
% [res, dat] = get_spatial_bootshuff(__, 'Name', Value) process with Name-Value pairs 
%
% INPUT
%
% 'pox'        - Numeric vector [N x 1] that specifies the animal's X coordinates.
%                Units are in cm. Required input.
%
% 'poy'        - Numeric vector [N x 1] that specifies the animal's Y coordinates.
%                Units are in cm. Required input.
%
% 'pot'        - Numeric vector [N x 1] that specifies the behavioral timestamps.
%                Units are in seconds. Required input.
%
% 'rmset'      - Structure containing settings for the rate map generation (see
%                rate_mapper)
%                Will be filled with defaults if empty
%
% 'sindx'      - (Name-Value) Numeric vector [M x 1] of indices mapping spikes to position samples.
%                Default is empty (must provide either sindx or spt).
%
% 'spt'        - (Name-Value) Numeric vector [M x 1] of spike timestamps.
%                Units are in seconds. Default is empty.
%
% 'iti'        - (Name-Value) Numeric vector [1 x 2] specifying the number of [bootstraps, shuffles].
%                Default value is [100 100].
%
% 'poh'        - (Name-Value) Numeric vector [N x 1] of head direction angles.
%                Units are in radians. Optional, but required for HD analysis.
%                Default is empty.
%
% 'srate'      - (Name-Value) Scalar, positive integer specifying the tracking sampling rate.
%                Units are in Hertz (Hz). Default value is 50.
%
% 'scores'     - (Name-Value) Cell array of strings specifying which analyses to run.
%                Options: 'all', 'spatial', 'grid', 'directional', 'speed',
%                'ahv', 'boundary'
%                Default is {'spatial'}.
%
% 'metrics'    - (Name-Value) Cell array of strings specifying which spatial info metrics to run.
%                Options: 'all', 'spatial_info', 'sparsity', 'mutual_info', etc. 
%                See get_spatial_info for all options.
%                Default is {'all'}.
%
% 'grid_type' - (Name-Value) String specifying the gridness algorithm to use.
%                Options: 'savelli', 'langston', 'wills', etc. 
%                See get_grid_score for all options.
%                Default is 'savelli'.
%
% 'trial'       - (Name-Value) Mx2 array of trial start times (column 1) and end
%                times (column 2), if this is provided, data within each trial
%                will be shuffled independently to other trial periods
%                Default is [pot(1), pot(end)]
%
% OUTPUT
%
% 'res'        - Table containing summary statistics for all calculated metrics.
%                Includes bootstrapped means, medians, 95% CIs, Z-Scores, and 
%                left/right-tailed p-values based on the null distribution.
%
% 'dat'        - Table containing the raw metric values for every bootstrap and shuffle iteration.
%                Column 1 ('IterationType') denotes 1=Bootstrap, 2=Shuffle.
%
% NOTES
% 1. The output table contains both a right-side p-value, and corresponding 95th
%    percentile. This should be used for metrics if you want to test that the
%    cell EXCEEDS the shuffles, for example spatial information content, where a
%    higher value than chance is good. The table also contains a left-sided
%    p-value and corresponding 5th percentile. This should be used for metrics
%    if you want to test that the cell UNDERPERFORMS chance, for example
%    sparsity, mutual information, shannon's entropy, where a low value is good.
% 
% 2. p-values are calculated as ((Nshuffles >= bootstrapped mean) + 1) / (Nshuffles + 1);
%    The +1 is to prevent a p-value of zero and is a standard correction for
%    this procedure, see: Phipson & Smyth (2010)
%
% 3. You must provide either 'sindx' or 'spt'. If 'spt' is provided, the function 
%    automatically calculates the nearest-neighbor indices to the position data.
%
% 4. Temporal shuffling uses a minimum circular offset of 20 seconds to ensure 
%    local spatial correlations are thoroughly broken in the null distribution.
%
% EXAMPLE
% 
%   % using get_synth_neurons
%   [trajectory,cells] = get_synth_neurons(1,0,0,"hpc_widths",[2 8]);
%   rmset.ppm = 100; % pixels per meter, can be scalar (same value usd for all sessions) or Nx1, were N = number of sesions
%   rmset.method = 'histogram';
%   rmset.binsize = 25;
%   rmset.ssigma = 30;
%   rmset.hd_type = 'histogram';
%   rmset.hd_bins = 60;
%   rmset.hd_boxcar = 3;
%   [res,dat] = get_spatial_bootshuff(trajectory.x(:),trajectory.y(:),trajectory.t(:),rmset,'spt',cells.PC(1).t,'poh',trajectory.hd(:),'scores',{'all'});
%   disp(res);
% 
% SEE ALSO get_spatial_info, get_grid_score, get_directional_info

% HISTORY
%
% version 1.0.0, Release 10/05/26 Initial release
% version 1.0.1, Release 11/05/26 Working release
% version 1.1.0, Release 07/07/26 Changed wrapping to account for periods not beginning at t=0
% version 1.2.0, Release 07/07/26 Atika fix: added ability to shuffle within trials
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
        pox {mustBeNumeric}
        poy {mustBeNumeric}
        pot {mustBeNumeric}
        rmset struct
        opts.sindx {mustBeNumeric} = []
        opts.spt {mustBeNumeric} = []
        opts.iti (1,2) double = [100 100]
        opts.poh {mustBeNumeric} = []
        opts.srate {mustBeNumeric} = 50
        opts.trial {mustBeNumeric} = []
        opts.epoly (:,2) double = []

        % what evaluation types to conduct
        opts.eval_type (1,:) string {mustBeMember(opts.eval_type, {'observed', 'shuffle', 'bootstrap'})} = ["observed", "shuffle", "bootstrap"]

        % what optional shuffles do we want to conduct
        opts.scores (1,:) string {mustBeMember(opts.scores, {'all','spatial','grid','directional','speed','boundary','ahv','anova'})} = ['spatial']

        % inputs passed to get_spatial_info
        opts.metrics (1,:) string {mustBeMember(opts.metrics, {'all','spatial_info','sparsity','mutual_info','entropy','kld','spatial_coherence','snr'})} = ['all']

        % inputs passed to get_grid_score
        opts.grid_type (1,1) string {mustBeMember(opts.grid_type, {'allen','wills','langston','soman','mixed','brandon','sargolini','krupic','savelli'})} = ['savelli']
    end
    if any(ismember(opts.scores,{'all'}))
        opts.scores = {'spatial','grid','directional','speed','boundary','ahv'};
    end

    % check spike index and generate if necessary
    if isempty(opts.sindx)
        if isempty(opts.spt)
            error('You must provide either a spike index or vector of spike times.')
        else
            opts.sindx = knnsearch(pot(:),opts.spt(:));
        end
    end
    opts.sindx = double(opts.sindx);

    % sort out spike data
    pos = [pox(:) poy(:)];
    spk = pos(opts.sindx,:);
    if isempty(spk)
        return
    end

    % head direction
    if ~isfield(opts,'poh') || isempty(opts.poh)
        if any(ismember(opts.scores,{'directional','boundary'}))
            warning('No head direction info provided, will skip directional & boundary analyses')
            nindx = ismember(opts.scores,{'directional','boundary'});
            opts.scores(nindx) = [];
        end
        poh = NaN(size(pox));
        sph = NaN(size(opts.sindx));
    else
        poh = opts.poh(:);
        sph = poh(opts.sindx);
    end

    % boundaries
    if any(ismember(opts.scores,{'boundary'}))
        if ~isfield(opts,'epoly') || isempty(opts.epoly)
            warning('No environment polygon provided, will skip boundary analyses')
            nindx = ismember(opts.scores,{'boundary'});
            opts.scores(nindx) = [];            
        end
    end

    % trial shuffle option
    if ~isfield(opts,'trial') || isempty(opts.trial)
        opts.trial = [pot(1), pot(end)];
    end
    n_trials = size(opts.trial, 1);

    % map trial start/end times to position indices (pot)
    trial_bounds = zeros(n_trials, 3); % [Start_Idx, End_Idx, Length]
    for k = 1:n_trials
        % find position indices that fall within this trial window
        idx_in_trial = find(pot >= opts.trial(k, 1) & pot <= opts.trial(k, 2));
        
        if ~isempty(idx_in_trial)
            trial_bounds(k, 1) = idx_in_trial(1);     % A_k (Start index)
            trial_bounds(k, 2) = idx_in_trial(end);   % B_k (End index)
            trial_bounds(k, 3) = numel(idx_in_trial); % L_k (Length in samples)
        end
    end
    valid_trials = trial_bounds(:, 3) > 0;
    trial_bounds = trial_bounds(valid_trials, :);
    n_trials = size(trial_bounds, 1);

    % preallocate
    speedlift = [];
    hd_speedlift = [];
    b_speedlift = struct;

    % precalculate shuffles
    offsets_matrix = zeros(n_trials, opts.iti(2));
    for k = 1:n_trials
        L_k = trial_bounds(k, 3);
        validOffsets = [-(L_k-1):-round(L_k/4), round(L_k/4):(L_k-1)]';
        
        % fallback if validOffsets is empty (for extremely short/empty trials)
        if isempty(validOffsets)
            validOffsets = 0;
        end
        
        % Draw random offsets for this specific trial across all shuffle iterations
        offsets_matrix(k, :) = validOffsets(randi(numel(validOffsets), opts.iti(2), 1));
    end
    % n = numel(pot); 

    % Setup ANOVA inputs
    if any(ismember(opts.scores, {'anova'}))
        [dwell_trials, bin_idx, total_bins] = setup_anova_bins(pox, poy, rmset, trial_bounds, opts.srate);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Bootstrap and shuffle
% This analysis has two steps:
% 1. Bootstrap (resample spikes with replacement)
%   The first is to calculate our 'observed' spatial
%   metrics, rather than just taking the initial value, to control for unstable or
%   unreliable activity, we will resample or bootstrap the firing rate map, for
%   example instead of calculating spatial information content once on the firing
%   rate map generated using the original spike train, we will resample the spike
%   train N times, generate N ratemaps, calculate N spatial information values and
%   then derive the mean of these. This provides a bootstrapped value which will
%   often not differ too much from the original, but which will decrease for
%   unstable cells.
% 2. Shuffle (shuffle spike train) 
%   Step 2 is to shuffle the spike train by randomly shifting it in time, this
%   preserves the 'local' spiking patterns and structure of the spike train, but
%   dissociates it from the animal's location. To deal with gaps in the data,
%   where recordings were paused or the data were maybe generated from
%   concatenated trials, we shuffle the spikes in the spike>position index
%   reference frame, this ensures spikes are not shuffled into missing time
%   periods.
%
% For simplicity, we will do both of these steps within the same loop, as
% resampling and shuffling differ only in the way the spikes are treated in each
% iteration. We will collect the data in dat and then calculate the result at
% the end.
    rng(999); % for reproducibility
    dat = table;

    % Build dynamic iteration queue based on user inputs
    % 1 = Bootstrap, 2 = Boot-Shuffle, 3 = Reg-Shuffle, 4 = Observed
    run_types = [];
    if ismember('bootstrap', opts.eval_type)
        run_types = [run_types, repmat(1, 1, opts.iti(1))]; 
        run_types = [run_types, repmat(2, 1, opts.iti(2))]; 
    end
    if ismember('shuffle', opts.eval_type)
        run_types = [run_types, repmat(3, 1, opts.iti(2))]; 
    end
    if ismember('observed', opts.eval_type) || ismember('shuffle', opts.eval_type)
        run_types = [run_types, 4]; % Always need 1 raw observed if shuffling or requested
    end

    % run the shuffles
    for iter = 1:length(run_types)
        current_type = run_types(iter);
        
        if current_type == 1 % Bootstrap (resample spikes with replacement)
            resample_idx = randi(length(opts.sindx), [length(opts.sindx), 1]); 
            sindx = opts.sindx(resample_idx); 
            spk_now = spk(resample_idx, :); 
            sph_now = sph(resample_idx);
            
        elseif current_type == 2 || current_type == 3 % Shuffle routines
            sindx = zeros(size(opts.sindx)); 
            
            for k = 1:n_trials
                A_k = trial_bounds(k, 1);
                B_k = trial_bounds(k, 2);
                L_k = trial_bounds(k, 3);
                offset_now = offsets_matrix(k, randi(opts.iti(2))); % Random offset from precalculated
                
                spk_mask = (opts.sindx >= A_k & opts.sindx <= B_k);
                if ~any(spk_mask)
                    continue 
                end
                
                sindx(spk_mask) = A_k + mod(opts.sindx(spk_mask) - A_k + offset_now, L_k);
            end
            
            valid_spikes = find(sindx > 0);
            
            % Apply resampling ONLY if it is a Boot-Shuffle (Type 2)
            if current_type == 2
                resample_idx = randi(length(valid_spikes), [length(valid_spikes), 1]);
                valid_spikes = valid_spikes(resample_idx);
            end
            
            sindx = sindx(valid_spikes);
            spk_now = pos(sindx, :);
            sph_now = poh(sindx, :);
            
        elseif current_type == 4 % Observed (raw data, no manipulation)
            sindx = opts.sindx;
            spk_now = spk;
            sph_now = sph;
        end

        if any( ismember(opts.scores,{'spatial','grid'}) )
            % bootstrapped firing rate map          
            [rmap_now,dmap_now,~,~,speedlift] = rate_mapper(pos.*10,spk_now.*10,rmset,speedlift);
        end
        
        % spatial metrics
        m = NaN(1,2);
        if any( ismember(opts.scores,{'spatial'}) )
            si = get_spatial_info(dmap_now,rmap_now,'metrics',opts.metrics);
            m = [si.skaggs_si_bits_per_sec, si.skaggs_si_bits_per_spike];                
        end

        % grid score
        g = NaN;
        if any( ismember(opts.scores,{'grid'}) )
            automap_boot = ndautoCORR(rmap_now,rmap_now,50);
            g = get_grid_score(automap_boot,rmset.binsize/10,'method',opts.grid_type); 
        end

        % head direction
        h = NaN(1,3);
        md = NaN(1,5);            
        if any( ismember(opts.scores,{'directional'}) )
            [hd_ratemap_now,hd_dwellmap_now,~,xi,~,hd_speedlift] = rate_mapper_hd(poh,sph_now,rmset,'speedlift',hd_speedlift);
            di = get_directional_info(hd_ratemap_now,hd_dwellmap_now,xi);
            h = [di.rayleigh_v, di.skaggs_si_bits_per_sec, di.skaggs_si_bits_per_spike];

            % multidirectionality
            [md, ~] = get_multidirectional_scores(hd_ratemap_now,[1:5]);
        end

        % angular head velocity
        ahv = NaN(1,4);
        if any( ismember(opts.scores,{'ahv'}) )
            [ahvscores, ~, ai] = get_ahv_info(pot, rad2deg(poh), sindx);
            ahv = [ahvscores,ai.skaggs_si_bits_per_sec, ai.skaggs_si_bits_per_spike];
        end

        % running speed
        sscore = NaN(1,3);
        if any( ismember(opts.scores,{'speed'}) )               
            [speedscore, ~, ~, si] = get_speed_info(pot, pox, poy, sindx);
            sscore = [speedscore, si.skaggs_si_bits_per_sec, si.skaggs_si_bits_per_spike];
        end

        % boundaries
        b = NaN(1,6);
        if any( ismember(opts.scores,{'boundary'}) )    
            [allo_stats, ego_stats, b_speedlift] = get_boundary_scores([pos poh], sindx, opts.epoly, opts.srate, 'speedlift',b_speedlift);
            b = [allo_stats.MRL, allo_stats.skaggs_si_bits_per_sec, allo_stats.skaggs_si_bits_per_spike,...
                ego_stats.MRL, ego_stats.skaggs_si_bits_per_sec, ego_stats.skaggs_si_bits_per_spike];
        end

        % Compute ANOVA F-statistic
        f_stat = NaN;
        if any(ismember(opts.scores, {'anova'}))
            f_stat = compute_anova_fstat(sindx, bin_idx, total_bins, n_trials, trial_bounds, dwell_trials);
        end

        % convert to a table
        i_t = array2table(current_type, "VariableNames", {'IterationType'});
        m_t = array2table(m,"VariableNames",{'XY_si_sec','XY_si_spike'});
        g_t = array2table(g,"VariableNames",{'grid_score'});
        h_t = array2table(h,"VariableNames",{'HD_rv','HD_si_sec','HD_si_spike'});
        md_t = array2table(md,"VariableNames",{'MD1','MD2','MD3','MD4','MD5'});
        a_t = array2table(ahv,"VariableNames",{'AHV_left','AHV_right','AHV_si_sec','AHV_si_spike'}); 
        s_t = array2table(sscore,"VariableNames",{'speed_score','speed_si_sec','speed_si_spike'}); 
        b_t = array2table(b,"VariableNames",{'BVC_MRL','BVC_si_sec','BVC_si_spike','EBC_MRL','EBC_si_sec','EBC_si_spike'}); 
        f_t = array2table(f_stat,"VariableNames",{'ANOVA_F'});   

        dat = [dat; i_t m_t, g_t, h_t, md_t, a_t, s_t, b_t, f_t];
    end

%%%%%%%%%%%%%%%% Get results
    metricNames = dat.Properties.VariableNames(2:end);
    
    bootMatrix = dat{dat{:, 1} == 1, 2:end};
    bootShufMatrix = dat{dat{:, 1} == 2, 2:end};
    regShufMatrix = dat{dat{:, 1} == 3, 2:end};
    obsMatrix = dat{dat{:, 1} == 4, 2:end};
    
    % Initialize outputs with NaNs in case certain evaluation types were skipped
    obsVals = NaN(1, size(obsMatrix, 2));
    if ~isempty(obsMatrix)
        obsVals = obsMatrix(1, :);
    end
    
    bootMean = mean(bootMatrix, 1, 'omitnan');
    bootMedian = median(bootMatrix, 1, 'omitnan'); 
    bootCI_low = prctile(bootMatrix, 2.5, 1); 
    bootCI_hi = prctile(bootMatrix, 97.5, 1); 
    
    % Observed values vs Regular Shuffle
    zScores_obs = NaN(size(obsVals));
    pValues_Right_obs = NaN(size(obsVals));
    pValues_Left_obs = NaN(size(obsVals));
    shuf95th = NaN(size(obsVals));
    shuf5th = NaN(size(obsVals));
    
    if ~isempty(regShufMatrix)
        numRegShuffles = size(regShufMatrix, 1);
        zScores_obs = (obsVals - mean(regShufMatrix, 1, 'omitnan')) ./ std(regShufMatrix, 0, 1, 'omitnan'); 
        
        pValues_Right_obs = (sum(regShufMatrix >= obsVals, 1) + 1) / (numRegShuffles + 1);
        pValues_Left_obs  = (sum(regShufMatrix <= obsVals, 1) + 1) / (numRegShuffles + 1);
        
        % Wipe p-values for metrics that were not computed
        nan_cols_obs = all(isnan(regShufMatrix), 1);
        pValues_Right_obs(nan_cols_obs) = NaN;
        pValues_Left_obs(nan_cols_obs) = NaN;
        
        shuf95th = prctile(regShufMatrix, 95, 1);
        shuf5th = prctile(regShufMatrix, 5, 1);
    end

    % Bootstrapped values vs Resampled Shuffle
    zScores = NaN(size(bootMean));
    pValues_Right = NaN(size(bootMean));
    pValues_Left = NaN(size(bootMean));
    
    if ~isempty(bootShufMatrix)
        numBootShuffles = size(bootShufMatrix, 1);
        zScores = (bootMean - mean(bootShufMatrix, 1, 'omitnan')) ./ std(bootShufMatrix, 0, 1, 'omitnan'); 
        
        pValues_Right = (sum(bootShufMatrix >= bootMean, 1) + 1) / (numBootShuffles + 1);
        pValues_Left  = (sum(bootShufMatrix <= bootMean, 1) + 1) / (numBootShuffles + 1);
        
        % Wipe p-values for metrics that were not computed
        nan_cols_boot = all(isnan(bootShufMatrix), 1);
        pValues_Right(nan_cols_boot) = NaN;
        pValues_Left(nan_cols_boot) = NaN;
    end

    % Compile everything into summary table
    res = table(obsVals', bootMean', bootMedian', bootCI_low', bootCI_hi', zScores', pValues_Right', shuf95th', pValues_Left', shuf5th', zScores_obs', pValues_Right_obs', pValues_Left_obs',...
                         'VariableNames', {'Observed', 'Boot_mean', 'Boot_median', 'Boot_CI_lower','Boot_CI_upper', 'z_score', 'p_value(right)', 'shuffle_95th_prctile','p_value(left)', 'shuffle_5th_prctile', 'z_score_obs', 'p_value(right)_obs','p_value(left)_obs'}, ...
                         'RowNames', metricNames);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%% setup_anova_bins
function [dwell_trials, bin_idx, total_bins] = setup_anova_bins(pox, poy, rmset, trial_bounds, srate)
    % Define spatial bin edges
    mid_x = mean(rmset.maplims([1 3]));
    mid_y = mean(rmset.maplims([2 4]));
    max_x = rmset.maplims([3]) - mid_x;
    max_y = rmset.maplims([4]) - mid_y;
    xvec = 0 : rmset.binsize : (max_x + rmset.binsize + rmset.padding); 
    yvec = 0 : rmset.binsize : (max_y + rmset.binsize + rmset.padding);    
    xvec = unique(sort([-xvec xvec], 'ascend')); 
    yvec = unique(sort([-yvec yvec], 'ascend'));
    
    % Discretize coordinates into 1D bin indices
    pox_c = pox - mid_x;
    poy_c = poy - mid_y;
    x_idx = discretize(pox_c, xvec);
    y_idx = discretize(poy_c, yvec);
    
    num_x_bins = length(xvec) - 1;
    num_y_bins = length(yvec) - 1;
    total_bins = num_x_bins * num_y_bins;
    
    valid_pos = ~isnan(x_idx) & ~isnan(y_idx);
    bin_idx = NaN(size(pox));
    bin_idx(valid_pos) = sub2ind([num_y_bins, num_x_bins], y_idx(valid_pos), x_idx(valid_pos));
    
    % Accumulate dwell times per trial
    n_trials = size(trial_bounds, 1);
    dwell_trials = zeros(n_trials, total_bins);
    for k = 1:n_trials
        idx_k = trial_bounds(k, 1):trial_bounds(k, 2);
        valid_k = bin_idx(idx_k);
        valid_k = valid_k(~isnan(valid_k));
        if ~isempty(valid_k)
            dwell_trials(k, :) = accumarray(valid_k(:), 1, [total_bins, 1])' * (1/srate);
        end
    end
    dwell_trials(dwell_trials < rmset.mindwell) = NaN;
end

%%%%%%%%%%%%%%%% compute_anova_fstat
function [f_stat] = compute_anova_fstat(sindx, bin_idx, total_bins, n_trials, trial_bounds, dwell_trials)
    % Aggregate spike counts per trial
    spike_trials = zeros(n_trials, total_bins);
    for k = 1:n_trials
        A_k = trial_bounds(k, 1);
        B_k = trial_bounds(k, 2);
        spk_mask = (sindx >= A_k & sindx <= B_k);
        
        if any(spk_mask)
            spk_bins = bin_idx(sindx(spk_mask));
            spk_bins = spk_bins(~isnan(spk_bins));
            if ~isempty(spk_bins)
                spike_trials(k, :) = accumarray(spk_bins(:), 1, [total_bins, 1])';
            end
        end
    end
    
    % Calculate sum of squares and F-statistic
    f_stat = NaN;
    rate_trials = spike_trials ./ dwell_trials;
    valid_rates = ~isnan(rate_trials);
    n_j = sum(valid_rates, 1); 
    
    valid_bins = n_j > 1; % Require at least two passes to compute variance
    
    if sum(valid_bins) > 1
        Y_j = sum(rate_trials(:, valid_bins), 1, 'omitnan') ./ n_j(valid_bins); 
        Y_tot = sum(rate_trials(:, valid_bins), 'all', 'omitnan') / sum(n_j(valid_bins)); 
        
        SSB = sum(n_j(valid_bins) .* (Y_j - Y_tot).^2);
        diff_sq = (rate_trials(:, valid_bins) - Y_j).^2;
        SSW = sum(diff_sq, 'all', 'omitnan');
        
        df_B = sum(valid_bins) - 1;
        df_W = sum(n_j(valid_bins)) - sum(valid_bins);
        
        if SSW > 0 && df_W > 0
            f_stat = (SSB / df_B) / (SSW / df_W);
        end
    end
end



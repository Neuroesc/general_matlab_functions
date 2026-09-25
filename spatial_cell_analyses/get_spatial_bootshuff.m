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

        % what optional shuffles do we want to conduct
        opts.scores (1,:) string {mustBeMember(opts.scores, {'all','spatial','grid','directional','speed','boundary','ahv'})} = ['spatial']

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
    n = numel(pot); 

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
% figure
% tiledlayout('flow')
    rng(999); % for reproducibility
    dat = table;
    for ii = 1:2
        for jj = 1:opts.iti(ii)
            if ii==1 % Bootstrap (resample spikes with replacement)
                sindx = randi(length(opts.sindx),[length(opts.sindx),1]); % new spike index
                spk_now = spk(sindx,:); % new spike x,y values
                sph_now = sph(sindx);
            elseif ii==2 % Shuffle (shuffle spike train) 
                sindx = zeros(size(opts.sindx)); % Pre-allocate shifted indices
                
                for k = 1:n_trials
                    A_k = trial_bounds(k, 1);
                    B_k = trial_bounds(k, 2);
                    L_k = trial_bounds(k, 3);
                    offset_now = offsets_matrix(k, jj);
                    
                    % find all spikes whose original position index falls inside this trial
                    spk_mask = (opts.sindx >= A_k & opts.sindx <= B_k);
                    
                    if ~any(spk_mask)
                        continue % Cell didn't fire during this trial, skip
                    end
                    
                    % apply interval-shifted modulo arithmetic
                    % s_new = A + mod(s - A + offset, L)
                    sindx(spk_mask) = A_k + mod(opts.sindx(spk_mask) - A_k + offset_now, L_k);
                end
                
                % Only keep spikes that successfully landed within a valid trial block
                valid_spikes = (sindx > 0);
                spk_now = pos(sindx(valid_spikes), :);
                sph_now = poh(sindx(valid_spikes), :);
            end
    
            if any( ismember(opts.scores,{'spatial','grid'}) )
                % bootstrapped firing rate map          
                [rmap_now,dmap_now,~,~,speedlift] = rate_mapper(pos.*10,spk_now.*10,rmset,speedlift);
            end
            
            % spatial metrics
            m = NaN(1,4);
            if any( ismember(opts.scores,{'spatial'}) )
                si = get_spatial_info(dmap_now,rmap_now,'metrics',opts.metrics);
                m = [si.skaggs_si_bits_per_sec, si.skaggs_si_bits_per_spike, si.kldivergence];                
            end
    
            % grid score
            g = NaN;
            if any( ismember(opts.scores,{'grid'}) )
                automap_boot = ndautoCORR(rmap_now,rmap_now,50);
                g = get_grid_score(automap_boot,rmset.binsize/10,'method',opts.grid_type); 
            end
    
            % head direction
            h = NaN(1,4);
            md = NaN(1,5);            
            if any( ismember(opts.scores,{'directional'}) )
                [hd_ratemap_now,hd_dwellmap_now,~,xi,~,hd_speedlift] = rate_mapper_hd(poh,sph_now,rmset,'speedlift',hd_speedlift);
                di = get_directional_info(hd_ratemap_now,hd_dwellmap_now,xi);
                h = [di.rayleigh_v, di.skaggs_si_bits_per_sec, di.skaggs_si_bits_per_spike, di.kldivergence];

                % multidirectionality
                [md, ~] = get_multidirectional_scores(hd_ratemap_now,[1:5]);
            end

            % angular head velocity
            ahv = NaN(1,5);
            if any( ismember(opts.scores,{'ahv'}) )
                [ahvscores, ~, ai] = get_ahv_info(pot, rad2deg(poh), sindx);
                ahv = [ahvscores,ai.skaggs_si_bits_per_sec, ai.skaggs_si_bits_per_spike, ai.kldivergence];
            end

            % running speed
            sscore = NaN(1,4);
            if any( ismember(opts.scores,{'speed'}) )               
                [speedscore, ~, ~, si] = get_speed_info(pot, pox, poy, sindx);
                sscore = [speedscore, si.skaggs_si_bits_per_sec, si.skaggs_si_bits_per_spike, si.kldivergence];
            end

            % boundaries
            b = NaN(1,8);
            if any( ismember(opts.scores,{'boundary'}) )    
                [allo_stats, ego_stats, b_speedlift] = get_boundary_scores([pos poh], sindx, opts.epoly, opts.srate, 'speedlift',b_speedlift);
                b = [allo_stats.MRL, allo_stats.skaggs_si_bits_per_sec, allo_stats.skaggs_si_bits_per_spike, allo_stats.kldivergence,...
                    ego_stats.MRL, ego_stats.skaggs_si_bits_per_sec, ego_stats.skaggs_si_bits_per_spike, ego_stats.kldivergence];
            end

            % convert to a table
            i_t = array2table(ii,"VariableNames",{'IterationType'});
            m_t = array2table(m,"VariableNames",{'XY_si_sec','XY_si_spike','XY_kld'});
            g_t = array2table(g,"VariableNames",{'grid_score'});
            h_t = array2table(h,"VariableNames",{'HD_rv','HD_si_sec','HD_si_spike','HD_kld'});
            md_t = array2table(md,"VariableNames",{'MD1','MD2','MD3','MD4','MD5'});
            a_t = array2table(ahv,"VariableNames",{'AHV_left','AHV_right','AHV_si_sec','AHV_si_spike','AHV_kld'}); 
            s_t = array2table(sscore,"VariableNames",{'speed_score','speed_si_sec','speed_si_spike','speed_kld'}); 
            b_t = array2table(b,"VariableNames",{'BVC_MRL','BVC_si_sec','BVC_si_spike','BVC_kld','EBC_MRL','EBC_si_sec','EBC_si_spike','EBC_kld'}); 
            
            dat = [dat; i_t m_t, g_t, h_t, md_t, a_t, s_t, b_t];

% if ii==2
%     if jj<30
% nexttile
% imagesc(rmap_now,'AlphaData',~isnan(rmap_now));
% daspect([1 1 1])
% colormap('turbo')
%     else
%         keyboard
%     end
% end


        end
    end

%%%%%%%%%%%%%%%% Get results
    % get metric names
    metricNames = dat.Properties.VariableNames(2:end);
    
    % split the data
    isBoot = dat{:, 1} == 1;
    isShuf = dat{:, 1} == 2;
    
    bootMatrix = dat{isBoot, 2:end};
    shufMatrix = dat{isShuf, 2:end};
    numShuffles = size(shufMatrix, 1);
    
    % calculate bootstrapped estimates
    bootMean = mean(bootMatrix, 1, 'omitnan');
    bootMedian = median(bootMatrix, 1, 'omitnan'); % Median is often safer for skewed metrics
    bootCI_low = prctile(bootMatrix, 2.5, 1); % Lower 95% Confidence bound
    bootCI_hi = prctile(bootMatrix, 97.5, 1); % Upper 95% Confidence bound
    
    % z-scores (parametric)
    shufMean = mean(shufMatrix, 1, 'omitnan');
    shufStd = std(shufMatrix, 0, 1, 'omitnan');
    zScores = (bootMean - shufMean) ./ shufStd;
    
    % significance (non-parametric)
    % P-value: what fraction of shuffles beat our bootstrapped mean?
    % Adding +1 to numerator and denominator is a standard correction 
    % to avoid p=0 exactly (pseudo-count).
    % RIGHT-TAILED (for spatial info, grid score, etc. - higher is better)
    pValues_Right = (sum(shufMatrix >= bootMean, 1) + 1) / (numShuffles + 1);
    shuf95th = prctile(shufMatrix, 95, 1);
    
    % LEFT-TAILED (For sparsity, entropy, etc. - lower is better)
    pValues_Left  = (sum(shufMatrix <= bootMean, 1) + 1) / (numShuffles + 1);
    shuf5th = prctile(shufMatrix, 5, 1);

    % Compile everything into summary table
    res = table(bootMean', bootMedian', bootCI_low', bootCI_hi',zScores', pValues_Right', shuf95th', pValues_Left', shuf5th',...
                         'VariableNames', {'Boot_mean', 'Boot_median', 'Boot_CI_lower','Boot_CI_upper', 'z_score', 'p_value(right)', 'shuffle_95th_prctile','p_value(left)', 'shuffle_5th_prctile'}, ...
                         'RowNames', metricNames);
    % dat.offset = zeros(size(dat,1),1);
    % dat.offset(isShuf) = offsets(:);






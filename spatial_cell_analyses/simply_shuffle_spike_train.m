function [spt2, spindx2] = simply_shuffle_spike_train(pot, spt, opts)
% simply_shuffle_spike_train shifts spike trains for shuffles
%
% Shifts the indices of the spikes relative to the position data, so spikes
% will always be shifted in increments of the position sampling interval.
% If trial boundaries are provided, spikes are circularly shuffled within 
% each trial block independently.
%
% USAGE
%
% spt2 = shift_spike_train(pot, spt) process with default settings
%
% [spt2, spindx2] = shift_spike_train(pot, spt) outputs the new spike index 
%
% [spt2, spindx2] = shift_spike_train(___, 'Name', Value) process with Name-Value pairs 
%
% INPUT
%
% 'pot'      - Numeric vector [N x 1]. Position data time stamps. 
%              Units are in seconds. Required.
% 
% 'spt'      - Numeric vector [M x 1]. Spike time stamps. 
%              Units are in seconds. Required.
%
% 'trial'    - (Name-Value) Numeric matrix [T x 2]. Start and end times for 
%              each trial. If provided, shuffling occurs independently within 
%              these blocks. Units are in seconds. Default is empty (whole session).
% 
% 'spindx'   - (Name-Value) Numeric vector [M x 1]. Nearest neighbor index 
%              in 'pot' for each spike in 'spt'.
%              Default value is dynamically computed via knnsearch(pot, spt).
% 
% 'minadd'   - (Name-Value) Numeric scalar. Minimum time the spike train must 
%              be shuffled. Units are in samples (e.g., 100 samples = 20s at 50Hz).
%              Default value is 100.
% 
% 'maxadd'   - (Name-Value) Numeric scalar. Maximum time the spike train can 
%              be shuffled. Units are in samples.
%              Default value is numel(pot) - 100.
%
% OUTPUT
%
% 'spt2'     - Numeric vector [K x 1]. New, shuffled spike times. 
%              Units are in seconds.
%
% 'spindx2'  - Numeric vector [K x 1]. For every spike in the shifted train, 
%              its new closest neighbor index in the position data.
%
% NOTES
% 1. Spikes are shuffled in increments of the position data sampling rate. 
%    This prevents spikes from landing in recording gaps.
%
% 2. Sub-sample spike timing is preserved by calculating the time difference 
%    between the original spike and its nearest frame, and re-applying that 
%    difference to the shifted frame.
%
% 3. If 'trial' boundaries are provided, any spikes occurring outside of these 
%    windows (e.g., during inter-trial intervals) are excluded from the output.
%
% EXAMPLE
% 
% % Shuffle spikes across the entire session
% spt2 = simply_shuffle_spike_train(pot, spt);
%
% % Shuffle spikes within independent trial blocks
% trial_times = [0, 120; 150, 300]; % Two trials
% spt2 = simply_shuffle_spike_train(pot, spt, 'trial', trial_times);
%
% SEE ALSO knnsearch, get_spatial_bootshuff
%
% HISTORY
% version 1.0.0, Release 06/11/19 Initial release
% version 2.0.0, Release 20/08/20 Changed to shifting indices
% version 3.0.0, Release 19/02/23 Updated inputs, renamed shift_spike_train
% version 4.0.0, Release 19/02/23 Updated comments
% version 5.0.0, Release 22/07/26 Renamed simply_shuffle_spike_train from shift_spike_train
% version 5.0.1, Release 22/07/26 Converted to arguments, fixed modulo, added trial blocks (Atika fix)
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
        pot {mustBeNumeric, mustBeVector}
        spt {mustBeNumeric, mustBeVector}
        opts.trial {mustBeNumeric} = []
        opts.spindx {mustBeNumeric} = []
        opts.minadd (1,1) double {mustBePositive, mustBeInteger} = 100
        opts.maxadd (1,1) double {mustBePositive, mustBeInteger} = numel(pot) - 100
    end
    
    % Get spike index if one was not provided
    if isempty(opts.spindx)
        opts.spindx = knnsearch(pot(:), spt(:));
    end

    % ensure indices are double precision to prevent uint saturation during negative shifts
    opts.spindx = double(opts.spindx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY   
    % handle empty trial matrix by assigning the entire session as 1 trial
    if isempty(opts.trial)
        opts.trial = [pot(1), pot(end)];
    end

    n_trials = size(opts.trial, 1);
    trial_bounds = zeros(n_trials, 3); % [Start_Idx, End_Idx, Length]

    % map trial start/end times to position indices (pot)
    for k = 1:n_trials
        idx_in_trial = find(pot >= opts.trial(k, 1) & pot <= opts.trial(k, 2));
        if ~isempty(idx_in_trial)
            trial_bounds(k, 1) = idx_in_trial(1);
            trial_bounds(k, 2) = idx_in_trial(end);
            trial_bounds(k, 3) = numel(idx_in_trial);
        end
    end

    % remove invalid/empty trials
    valid_trials = trial_bounds(:, 3) > 0;
    trial_bounds = trial_bounds(valid_trials, :);
    n_trials = size(trial_bounds, 1);

    % pre-allocate shifted indices
    spindx2_full = zeros(size(opts.spindx)); 
    
    % block-wise modular shifting
    for k = 1:n_trials
        A_k = trial_bounds(k, 1);
        B_k = trial_bounds(k, 2);
        L_k = trial_bounds(k, 3);
        
        % dynamically adjust max/min shift based on this specific trial's length
        max_shift = min(opts.maxadd, L_k - 1);
        
        if L_k <= (2 * opts.minadd)
            % trial is too short for the requested minadd, fall back to 1/4 trial length
            safe_min = max(1, round(L_k / 4));
            validOffsets = [-(max_shift):-safe_min, safe_min:(max_shift)]';
        else
            validOffsets = [-(max_shift):-opts.minadd, opts.minadd:(max_shift)]';
        end
        
        % fallback if validOffsets is completely empty
        if isempty(validOffsets)
            validOffsets = 0;
        end
        
        % randomly pick exactly ONE shift value for this trial
        offset_now = validOffsets(randi(numel(validOffsets)));
        
        % find all spikes whose original position index falls inside this trial
        spk_mask = (opts.spindx >= A_k & opts.spindx <= B_k);
        
        if any(spk_mask)
            % apply interval-shifted modulo arithmetic: s_new = A + mod(s - A + offset, L)
            spindx2_full(spk_mask) = A_k + mod(opts.spindx(spk_mask) - A_k + offset_now, L_k);
        end
    end
    
    % extract only the valid spikes (filters out ITI spikes)
    valid_mask = spindx2_full > 0;
    spindx2 = spindx2_full(valid_mask);
    
    % preserve exact sub-frame spike timing
    % time of original spike minus time of original nearest frame
    spike_diff = spt(valid_mask) - pot(opts.spindx(valid_mask));
    
    % apply to new frame time
    spt2 = pot(spindx2) + spike_diff;

end
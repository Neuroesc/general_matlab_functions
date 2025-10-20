function [spindx2,spt2] = shuffle_spike_train(pot,spt,varargin)
%shiftSTRAIN  shifts spike trains for shuffles
% Shifts the indices of the spikes relative to the position data, so spikes
% will always be shifted in increments of the positions sampling interval.
%
%   spt2 = shift_spike_train(pot,spt) process with default settings
%   randomly shuffles vector of spike times spt by between 20s and max(spt) minus 20s
%
%   [spt2,spindx2] = shift_spike_train(pot,spt) outputs a new spike index (see below)
%
%   [spt2,spindx2] = shift_spike_train(___,Name,Value,...) Name-Value pairs can be used 
%   to control aspects of the shuffle, see below
%
%   parameters include:
%
%   'pot'           -   Numeric vector, required, position data time stamps, units are in seconds.
% 
%   'spt'           -   Numeric vector, required, spike time stamps, units are in seconds.
% 
%   'spindx'        -   Numeric vector, nearest neighbor in pot for each spike time in spt
%                       Default value is knnsearch(pot,spt)
% 
%   'minadd'        -   Numeric scalar, minimum time the spike train must be shuffle, units are 
%                       in samples not seconds
%                       Default value is 100 (20s at 50Hz sampling rate)
% 
%   'maxadd'        -   Numeric scalar, maximum time the spike train can be shuffle, units are 
%                       in samples not seconds
%                       Default value is numel(pot)-100 (20s at 50Hz sampling rate)
%
%   'potn'          -   Numeric scalar, total number of position data time stamps.
%                       Must be provided if pot is not provided.
%
%   outputs include:
%
%   'spt2'          -   [Nx1] Numeric vector, same size as spt, new spike times, units are in seconds
%
%   'spindx2'       -   [Nx1] Numeric vector, same size as spt, for every spike in the shifted spike train
%                       the new closest neighbor in the position data
%
%   Class Support
%   -------------
%   The input matrices pot and spt must be a real, non-sparse matrix of
%   the following classes: uint8, int8, uint16, int16, uint32, int32,
%   single or double.
%
%   Notes
%   -----
%   1. Spikes are shuffled in increments of position data sampling, this is much faster than using an arbitrary time value
%      and this means spikes are only shifted to possible position samples, which works even when there are gaps in the 
%      position time vector, which can happen when recordings are stitched together and we want to shuffle across them. As
%      most analyses, such as firing rate maps, only care about which position data sample a spike corresponds to, in most
%      (spatial) analyses this shuffle is sufficient.
%
%   Example
%   ---------
% 
%   % shuffle spikes using default values
%   pot = (1:50000)'.*(1/50); % dummy position data times
%   spt = sort(1000.*rand(100,1)); % dummy spike times
%   spindx = knnsearch(pot,spt);
%   spt2 = shift_spike_train(pot,spt);
%
%   % shuffle spikes with all the default values specified
%   spt2 = shift_spike_train(pot,spt,'spindx',spindx,'minadd',100,'maxadd',numel(pot)-100);
%
%   See also knnsearch

% HISTORY
%
% version 1.0.0, Release 06/11/19 Initial release
% version 2.0.0, Release 20/08/20 Changed to shifting indices
% version 3.0.0, Release 19/02/23 Updated inputs, renamed shift_spike_train
% version 4.0.0, Release 23/07/25 Updated comments
% version 4.1.0, Release 23/07/25 Renamed shuffle_spike_train to avoid conflicts
% version 4.2.0, Release 23/07/25 Improved inputs and removed requirement for spike/position times
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    p = inputParser;
    addParameter(p,'pot',@(x) ~isempty(x) && ~all(isnan(x(:))) );      
    addParameter(p,'spt',@(x) ~isempty(x) && ~all(isnan(x(:))) );  
    addParameter(p,'spindx',[],@(x) ~isempty(x) && ~all(isnan(x(:))) )
    addParameter(p,'minadd',100,@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'maxadd',[],@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'potn',[],@(x) isnumeric(x) && isscalar(x));     
    parse(p,pot,spt,varargin{:});
    config = p.Results;

    % get spike index if one was not provided
    if isempty(config.spindx) || all(isnan(config.spindx))
        if isempty(config.pot) || all(isnan(config.pot)) || isempty(config.spt) || all(isnan(config.spt))
            error('if a spike index (spindx) is not provided spike and position times (pot & spt) must be provided')
        else
            config.spindx = knnsearch(pot,spt);
        end
    end
    config.spindx = double(config.spindx);

    % calculate potn if it was not provided
    if isempty(config.potn) || all(isnan(config.potn))
        if isempty(config.pot) || all(isnan(config.pot))
            error('if total position samples (potn) is not provided position times (pot) must be provided')
        else
            config.potn = numel(pot);
        end
    end

    % calculate maxadd if it was not provided
    if isempty(config.maxadd) || all(isnan(config.maxadd))
        if isempty(config.potn) || all(isnan(config.potn))
            error('total position samples (potn) is required to calculate maximum time offset (maxadd)')
        else
            config.maxadd = config.potn-100;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % shift spike indices
    shift_amount = randi([config.minadd,config.maxadd],1) .* (2*(rand(1)>.5)-1);
    spindx2 = config.spindx + shift_amount;

    % wrap indices so they do not exceed session length
    spindx2(spindx2<1) = spindx2(spindx2<1)+config.potn;
    spindx2(spindx2>config.potn) = spindx2(spindx2>config.potn)-config.potn;

    % calculate new spike times if required
    if nargout>1
        spt2 = [];
        if ~isempty(config.pot) && ~all(isnan(config.pot)) && ~isempty(config.spt) && ~all(isnan(config.spt))
            spike_diff = config.pot(config.spindx)-config.spt; % difference between each spike and its nearest position time neighbour
            spt2 = config.pot(spindx2)-spike_diff;
        end
    end






























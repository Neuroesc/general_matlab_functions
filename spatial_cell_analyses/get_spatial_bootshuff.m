function b = get_spatial_bootshuff(pox,poy,pot,rmset,opts)
% name short description
% longer description
%
% USAGE
%
% out = name(in) process with default settings
%
% out = name(in,optional) process using optional argument 1
%
% out = rate_mapper(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'in' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% OUTPUT
%
% 'out' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% NOTES
% 1. 
%
% 2. 
%
% EXAMPLE
% 
% SEE ALSO NAME, NAME

% HISTORY
%
% version 1.0.0, Release 00/00/26 Initial release
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
        opts.iti double = 100
        opts.poh {mustBeNumeric} = []

        % what optional shuffles do we want to conduct
        opts.scores (1,:) string {mustBeMember(opts.scores, ["all","spatial","grid","directional"])} = ["all"]

        % inputs passed to get_spatial_info
        opts.metrics (1,:) string {mustBeMember(opts.metrics, ["all","spatial_info","sparsity","mutual_info","entropy","kld","spatial_coherence","snr"])} = ["all"]

        % inputs passed to get_grid_score
        opts.grid_score (1,1) string {mustBeMember(opts.grid_score, ['allen','wills','langston','soman','mixed','brandon','sargolini','krupic','savelli'])} = ["savelli"]
    end

    % check spike index and generate if necessary
    if isempty(opts.sindx)
        if isempty(opts.spt)
            error('You must provide either a spike index or vector of spike times.')
        else
            opts.sindx = knnsearch(pot,opts.spt);
        end
    end

    % sort out spike data
    pos = [pox(:) poy(:)];
    spk = pos(opts.sindx,:);
    if isempty(spk)
        return
    end

    % head direction
    if any(ismember({"directional","all"},opts.scores))
        poh = opts.poh(:);
        sph = poh(opts.sindx);
    end

    % preallocate
    speedlift = [];
    hd_speedlift = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Bootstrap (resample spikes with replacement)
% This analysis has two steps, the first is to calculate our 'observed' spatial
% metrics, rather than just taking the initial value, to control for unstable or
% unreliable activity, we will resample or bootstrap the firing rate map, for
% example instead of calculating spatial information content once on the firing
% rate map generated using the original spike train, we will resample the spike
% train N times, generate N ratemaps, calculate N spatial information values and
% then derive the mean of these. This provides a bootstrapped value which will
% often not differ too much from the original, but which will decrease for
% unstable cells.
    rng(999); % for reproducibility
    bootstrapped_data = table;
    for ii = 1:iti
        idx = randi(length(spindx),[length(spindx),1]); % new spike index
        spk_now = spk(idx,:); % new spike x,y values

        if any(ismember({"spatial","grid","all"},opts.scores))
            % bootstrapped firing rate map          
            [ratemap_boot,~,~,~,speedlift] = rate_mapper(pos,spk_now,rmset,speedlift);
        end
        
        % spatial metrics
        m = struct;
        if any(ismember({"spatial","all"},opts.scores))
            m = get_spatial_info(dmap,rmap,'metrics',opts.metrics);
        end

        % grid score
        g = struct;
        if any(ismember({"grid","all"},opts.scores))
            automap_boot = ndautoCORR(ratemap_boot,ratemap_boot,50);
            g = get_grid_score(automap_boot,rmset.binsize/10,'method',opts.grid_score); 
        end

        % head direction
        h = struct;
        if any(ismember({"directional","all"},opts.scores))
            sph_now = sph(idx);
            [hd_ratemap_boot,~,~,xi,~,hd_speedlift] = rate_mapper_hd(poh,sph_now,rmset,'dwellmap',hd_speedlift);
            h = get_directional_info(hd_ratemap_boot,xi);
        end

        % convert to a table
        bootstrapped_data = [bootstrapped_data; struct2table([m; g; h])];
    
    end  












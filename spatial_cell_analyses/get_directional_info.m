function h = get_directional_info(rmap,xi,opts)
% get_directional_info calculate various head direction measures
% Function calculates the most commonly used head direction measures and
% statistics, using many functions from the circular statistics toolbox
%
% USAGE
%
% m = get_directional_info(dmap,rmap,xi) process with default settings
%
% m = get_directional_info(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'rmap' - 1xM numeric matrix, required, firing rate map
%       units are in Hz (spikes per second).
%
% 'xi' - 1xM, bin locations, required
%       units are in radians
%
% 'metrics' - Cell array of strings, optional
%       Choices include: 
%           {'all'} - compute all spatial metrics
%           {'rayleigh'} - compute Rayleigh vector length
%           {'stats'} - compute PFD (max), mean, standard deviation
%       default is {'all'}
%
% OUTPUT
%
% 'm' - Structure, fields will depend on the outputs requested:
%       m.rv = Rayleigh vector length
%       m.pfd = Direction of maximum or preferred firing direction
%       m.mean = Circular mean
%       m.stdev = Circular standard deviation
%
% NOTES
% 1. See get_spatial_info, this can be used to calculate spatial information
% content, Kullback–Leibler divergence and other information theoretic values
% and can work with head direction dwell/rate maps
%
% EXAMPLE
%
% 
% SEE ALSO get_spatial_info

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
        rmap double
        xi double
        opts.metrics (1,:) string {mustBeMember(opts.metrics, ["all", "rayleigh", "stats"])} = ["all"]
    end
    h = struct;

    % prep
    rmap = rmap(:);
    xi = xi(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Circular statistics
    if ismember("rayleigh",opts.metrics) || ismember("all",opts.metrics)          
        h.rv = circ_r(xi,rmap); % rayleigh vector length
    end

%%%%%%%%%%%%%%%% Stats
    if ismember("stats",opts.metrics) || ismember("all",opts.metrics)          
        h.pfd = rad2deg( xi(hd3n == max(rmap)) ); % preferred angle (location of max frate)
        h.mean = rad2deg( circ_mean(xi,rmap) ); % mean angle
        h.stdev = rad2deg( circ_std(xi,rmap) ); % std deviation angle
    end

end













































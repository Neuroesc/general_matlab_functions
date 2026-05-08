function [ratemap,dwellmap,spikemap,xi,rmset,speedlift] = rate_mapper_hd(pos,spk,rmset,opts)
% rate_mapper_hd map directional spike and position data in 2D
%
% USAGE
%
% [ratemap,dwellmap,spikemap,rmset] = rate_mapper_hd(pos,spk,rmset) process with default settings
%
% [ratemap,dwellmap,spikemap,rmset] = rate_mapper_hd(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'pos' - Nx1, numeric matrix, the position head direction values
%       units are in radians!
%
% 'spk' - Nx1, numeric matrix, the spike head direction values data
%       units are in radians!
%
% 'rmset' - Struct, optional, additional settings, which can include:
%       rmset.hd_bins = Scalar, number of bins to use when binning data, default is 60
%       rmset.srate = Scalar, sampling rate of the position data, in Hz, default is 50
%       rmset.hd_sigma = Scalar, for density method, smoothing strength, default is 0.04
%       rmset.hd_boxcar = Integer, for histogram method, smoothing strength, default is 3
%       rmset.mindwell = Scalar, bins with a dwell time value less than this will be set to NaN, default is 0.01
%
% OUTPUT
%
% 'ratemap' - 1xM, numeric vector, firing rate map
%       units are in Hz.
%
% 'dwellmap' - 1xM, numeric vector, dwell time map
%       units are in seconds.
%
% 'spikemap' - 1xM, numeric vector, spike map
%       units are in spikes
%
% 'xi' - 1xM, bin locations
%       units are in radians
%
% 'rmset' - Struct, see inputs
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
% version 1.0.0, Release 08/05/26 Initial release, created from mapHD
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
        pos {mustBeNumeric} % radians
        spk {mustBeNumeric} % radians
        rmset struct = []
        opts.speedlift {mustBeNumeric} = []
        opts.type (1,1) string {mustBeMember(opts.type, ["histogram", "density"])} = ["histogram"]
    end

    if max(abs(pos),[],"all",'omitmissing')>2*pi
        warning('Head direction angles may not be in radians, please check')
    end

    pos = wrapTo2Pi(pos(:));
    spk = wrapTo2Pi(spk(:));

    % defaults in rmset
    defset = struct;
    defset.hd_sigma = 0.04;
    defset.hd_boxcar = 3;    
    defset.hd_bins = 60;
    defset.mindwell = 0.01;
    defset.srate = 50;

    % Fill in missing inputs in rmset using defset
    f1 = fieldnames(defset);
    if ~exist('rmset','var')
        rmset = struct;
    end
    f2 = fieldnames(rmset);
    for i = 1:size(f1,1)
        if ~ismember(f1{i},f2) % if this defset field does not exist in rmset
            rmset.(f1{i}) = defset.(f1{i}); % add it to rmset
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY   
    % Head direction analysis
    xi = linspace(0,2*pi,rmset.hd_bins)'; % angles for binning

    switch opts.type
        case {'density'}
            % Create a firing rate by head direction map using a kernel smoothed density estimate approach
            % Use a head direction dwell map if one was provided, or create it if not
            if isempty(opts.speedlift) || all(isnan(opts.speedlift(:)))
                [hd1] = circ_ksdensity(pos,xi,[],rmset.hd_sigma); % the session head direction       
                dwellmap = hd1 .* (1/50); % convert session HD to time 
            else
                dwellmap = opts.speedlift;
            end
            speedlift = dwellmap;

            % Calculate the kernel smoothed spike map and then the firing rate map
            [spikemap] = circ_ksdensity(spk,xi,[],rmset.hd_sigma); % the cell's head direction
            ratemap = spikemap ./ dwellmap; % calculate HD firing rate
            
        case {'histogram'}
            % Create a firing rate by head direction map using a straight histogram approach
            % Use a head direction dwell map if one was provided, or create it if not
            fh = fspecial('average',[1 rmset.hd_boxcar]); % boxcar filter                
            if isempty(opts.speedlift) || all(isnan(opts.speedlift(:)))
                hd1 = histcounts(pos,xi); % the session head direction   
                hd1 = imfilter(hd1,fh,'circular','same');         
                hc = movmean(xi,2,'EndPoints','discard');
                dwellmap = interp1([hc-2*pi; hc; hc+2*pi],[hd1(:);hd1(:);hd1(:)],xi,'linear'); % triplicate and interpolate
                dwellmap = dwellmap .* (1/50); % convert session HD to time
            else
                dwellmap = opts.speedlift;
            end            
            speedlift = dwellmap;

            % Calculate the binned spike map and then the firing rate map            
            spikemap = histcounts(spk,xi); % the cell's head direction
            spikemap = imfilter(spikemap,fh,'circular','same'); % smoothe the data  
            hc = movmean(xi,2,'EndPoints','discard');
            spikemap = interp1([hc-2*pi; hc; hc+2*pi],[spikemap(:);spikemap(:);spikemap(:)],xi,'linear'); % triplicate and interpolate
            ratemap = spikemap ./ dwellmap; % calculate HD firing rate
    end

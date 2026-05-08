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
        opts.grid_score {mustBeNumeric} = 0        
        opts.rayleigh_vector {mustBeNumeric} = 0

        % inputs passed to get_spatial_info
        opts.spatial_info {mustBeNumeric} = 0
        opts.sparsity {mustBeNumeric} = 0
        opts.mutual_info {mustBeNumeric} = 0
        opts.entropy {mustBeNumeric} = 0
        opts.kld {mustBeNumeric} = 0
        opts.spatial_coherence {mustBeNumeric} = 0
        opts.snr {mustBeNumeric} = 0
        opts.all {mustBeNumeric} = 1
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
    speedlift = [];
    hd_dwellmap = [];


    rng(999); % for reproducibility
    for ii = 1:iti
        idx = randi(length(spindx),[length(spindx),1]); % new spike index
        spk_now = spk(idx,:);

        % bootstrapped firing rate map          
        [ratemap_boot,~,~,~,speedlift] = rate_mapper(pos,spk_now,rmset,speedlift);
        
        % spatial metrics
        m = get_spatial_info(dmap,rmap,'all',opts.all,'spatial_info',opts.spatial_info,'sparsity',opts.sparsity,);


        % % Skaggs spatial information content (bits per second)
        % pr = dwellmap_boot ./ sum(dwellmap_boot(:),'omitnan'); % dwell time probability
        % ro = sum(ratemap_boot(:) .* pr(:),'omitnan'); % overall firing rate
        % si = sum(pr(:) .* (ratemap_boot(:)./ro) .* log2(ratemap_boot(:)./ro),'omitnan'); 

        % autocorrelation & grid score
        automap_boot = ndautoCORR(ratemap_boot,ratemap_boot,50);
        [g,~] = get_grid_score(automap_boot,rmset.binsize/10,'method','savelli'); 

        % % rayleigh vector
        % sph_now = sph(idx);
        % [~,hd_dwellmap,~,~,r,~,~,~] = mapHD(hd_dwellmap,poh(:),sph_now(:),rmset);
        
        % accumulate
        % shuff_info.si_boot(ii) = si;
        shuff_info.g_boot(ii) = g;
        % shuff_info.r_boot(ii) = r;        
    end  













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JointEntropy sub function
function [hd_rmap,hd_dmap] = hd_maps(hd_dmap,poh,sph,rmset) % in radians
    % Head direction analysis
    ai = linspace(0,2*pi,rmset.hd_bins)'; % angles for binning

    % Create a firing rate by head direction map using a straight histogram approach
    % Use a head direction dwell map if one was provided, or create it if not
    fh = fspecial('average',[1 config.hd_boxcar]); % boxcar filter                
    if ~exist('hd_dmap','var') || isempty(hd_dmap) || all(isnan(hd_dmap(:)))
        hd1 = histcounts(poh,ai); % the session head direction   
        hd1 = imfilter(hd1,fh,'circular','same');         
        hc = movmean(ai,2,'EndPoints','discard');
        hd_dmap = interp1([hc-2*pi; hc; hc+2*pi],[hd1(:);hd1(:);hd1(:)],ai,'linear'); % triplicate and interpolate
        hd_dmap = hd_dmap .* (1/50); % convert session HD to time
    end              
    
    % Calculate the binned spike map and then the firing rate map            
    spikemap = histcounts(sph,ai); % the cell's head direction
    spikemap = imfilter(spikemap,fh,'circular','same'); % smoothe the data  
    hc = movmean(ai,2,'EndPoints','discard');
    spikemap = interp1([hc-2*pi; hc; hc+2*pi],[spikemap(:);spikemap(:);spikemap(:)],ai,'linear'); % triplicate and interpolate
    hd_rmap = spikemap ./ hd_dmap; % calculate HD firing rate

end





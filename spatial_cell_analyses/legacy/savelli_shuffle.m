function [shuff_info] = savelli_shuffle(pox,poy,pot,poh,spt,spindx,rmset,iti)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% savelli_shuffle  compute grid score or spatial information with shuffle 
% Analysis from: https://elifesciences.org/articles/21354, Savelli et al. (2017)
%
% "Visual inspection of rate maps suggested to us that the exclusive use of a single 
% gridness score threshold, however determined, could not keep the rate of both false 
% positives and false negatives at a satisfactory level in our dataset and for our 
% study’s goals. Our analyses were particularly sensitive to the accuracy of the 
% estimation of grid parameters, but we did not find the gridness score to provide 
% a reliable measure of how ‘clean’ the grid was."
% This analysis includes a bootstrapping procedure where grid score is
% calculated on 100 bootstrapped firing rate maps, the the grid score then 
% has to be ≥0.1 for at least 95 out of the 100 bootstrapped rate maps.
%
% USAGE:
%       [out] = template(in) process with default settings
% 
%       [out] = template(in,optional1) process using optional argument 1
% 
%       [out] = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'param1'          -   (default = X) Scalar value, parameter to do something
% 
%       'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%       in    - input as a vector
%
% OUTPUT:
%       out   - output as a vector
%
% EXAMPLES:
%       % run function using default values
%       out = template(in,varargin)
%
% See also: GIT_audit

% HISTORY:
% version 1.0.0, Release 16/02/23 Code conception
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    % default settings for creating firing rate maps
    defset              = struct;
    defset.method       = 'histogram';
    defset.binsize      = 20;
    defset.ssigma       = 80;
    defset.ash          = 16;
    defset.maplims      = [min(pox(:)) min(poy(:)) max(pox(:)) max(poy(:))];     
    defset.padding      = 0; % in mm
    defset.mindwell     = 0; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    defset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    defset.maxdist      = 640; % (mm, default 50) used by kadaptive, adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    defset.srate        = 50; % (default 50) sampling rate of data in Hz, used to calculate time
    defset.steps        = 32; % the number of convolution size steps to use for kadaptive    
    defset.kern         = 'biweight'; % kernel
    defset.smethod      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing
    defset.bmethod      = 0; % binning method, 0 = use rmset.binsize and create a grid with binsize x binsize pixels, N = create a grid that is NxN and ends at the data limits

    %% check that all parameters are included in rmset 
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

    % iterations (number of shuffles)
    if ~exist('iti','var')
        iti = 100;
    end     

    spx = pox(spindx);
    spy = poy(spindx);
    sph = poh(spindx);    
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Bootstrap (resample spikes with replacement)  
    shuff_info = array2table(NaN(iti,6,'single'),'VariableNames',{'si_boot','g_boot','r_boot','si_shuff','g_shuff','r_shuff'});
    if isempty(spx)
        return
    end
    speedlift = [];
    hd_dwellmap = [];
    pos = [pox(:) poy(:)];
    spk = [spx(:) spy(:)];
    poh = rad2deg(poh(:));
    sph = rad2deg(sph(:));
    rng(999); % for reproducibility
    for ii = 1:iti
        idx = randi(length(spindx),[length(spindx),1]); % new spike index
        spk_now = spk(idx,:);

        % bootstrapped firing rate map          
        [ratemap_boot,~,~,~,speedlift] = rate_mapper(pos,spk_now,rmset,speedlift);
        
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
    
%% >>>>>>>>>> Shuffle (shuffle spike train) 
    rng(999); % for reproducibility
    n = numel(pot); 
    minShift = rmset.srate*20; % exclude offsets less than 20s
    validOffsets = [-n:-minShift, minShift:n]; % Build the allowed offset set
    offsets = validOffsets(randi(numel(validOffsets),iti,1)); % Randomly pick one 

    for ii = 1:iti
        % [~,spindx2] = shift_spike_train(pot,spt,'spindx',spindx);
        spindx2 = mod(spindx-1 + offsets(ii), n) + 1;
        spk_now = pos(spindx2,:);

        % bootstrapped firing rate map   
        [ratemap_shuff,dwellmap_shuff,~,~,speedlift] = rate_mapper(pos,spk_now,rmset,speedlift);

        % Skaggs spatial information content (bits per second)
        pr = dwellmap_shuff ./ sum(dwellmap_shuff(:),'omitnan'); % dwell time probability
        ro = sum(ratemap_shuff(:) .* pr(:),'omitnan'); % overall firing rate
        si = sum(pr(:) .* (ratemap_shuff(:)./ro) .* log2(ratemap_shuff(:)./ro),'omitnan'); 

        % autocorrelation & grid score
        automap_shuff = ndautoCORR(ratemap_shuff,ratemap_shuff,50);
        [g,~] = get_grid_score(automap_shuff,rmset.binsize/10,'method','savelli'); 

        % rayleigh vector
        sph_now = poh(spindx2,:);
        [~,hd_dwellmap,~,~,r,~,~,~] = mapHD(hd_dwellmap,poh(:),sph_now(:),rmset);
        
        % accumulate
        shuff_info.si_shuff(ii) = si;
        shuff_info.g_shuff(ii) = g;
        shuff_info.r_shuff(ii) = r;                
    end   












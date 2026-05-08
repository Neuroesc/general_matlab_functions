function [ai,dwellmap,spikemap,ratemap,r,mx,mn,sd] = mapHD(dwellmap,ppoh,psph,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         in - input 1
%         in2 - input 2
%
% OUTPUT:
%    p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 11/04/19 Initial release, created to contain this code instead of klustest
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
    narginchk(4,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    ppoh = ppoh(:);
    psph = psph(:);
    if any(ppoh(:)<0)
        ppoh = wrapTo360(ppoh);
        psph = wrapTo360(psph);        
    end
    
    % Head direction analysis
    ai = linspace(0,2*pi,config.hd_bins)'; % angles for binning
    hd_c = deg2rad(psph); % cell head direction in radians

    switch config.hd_type
        case {'density'}
            % Create a firing rate by head direction map using a kernel smoothed density estimate approach
            % Use a head direction dwell map if one was provided, or create it if not
            if ~exist('dwellmap','var') || isempty(dwellmap) || all(isnan(dwellmap(:)))
                hd_s = deg2rad(ppoh); % session head direction in radians
                hd_s = wrapTo360(hd_s);                
                [hd1] = circ_ksdensity(hd_s,ai,[],config.hd_sigma); % the session head direction       
                dwellmap = hd1 .* (1/50); % convert session HD to time 
            end

            % Calculate the kernel smoothed spike map and then the firing rate map
            [spikemap] = circ_ksdensity(hd_c,ai,[],config.hd_sigma); % the cell's head direction
            ratemap = spikemap ./ dwellmap; % calculate HD firing rate

            % head direction analyses
            hd3n = ratemap ./ max(ratemap); % normalise cell hd
            hd3n = hd3n(:);            
            r = circ_r(ai,hd3n(:)); % rayleigh vector length
            mx = rad2deg(ai(hd3n == max(hd3n))); % preferred angle (location of max frate)
            mn = rad2deg(circ_mean(ai,hd3n)); % mean angle
            sd = rad2deg(circ_std(ai,hd3n)); % std deviation angle
            
        case {'histogram'}
            % Create a firing rate by head direction map using a straight histogram approach
            % Use a head direction dwell map if one was provided, or create it if not
            fh = fspecial('average',[1 config.hd_boxcar]); % boxcar filter                
            if ~exist('dwellmap','var') || isempty(dwellmap) || all(isnan(dwellmap(:)))
                hd1 = histcounts(deg2rad(ppoh),ai); % the session head direction   
                hd1 = imfilter(hd1,fh,'circular','same');         
                hc = movmean(ai,2,'EndPoints','discard');
                dwellmap = interp1([hc-2*pi; hc; hc+2*pi],[hd1(:);hd1(:);hd1(:)],ai,'linear'); % triplicate and interpolate
                dwellmap = dwellmap .* (1/50); % convert session HD to time
            end              

            % Calculate the binned spike map and then the firing rate map            
            spikemap = histcounts(deg2rad(psph),ai); % the cell's head direction
            spikemap = imfilter(spikemap,fh,'circular','same'); % smoothe the data  
            hc = movmean(ai,2,'EndPoints','discard');
            spikemap = interp1([hc-2*pi; hc; hc+2*pi],[spikemap(:);spikemap(:);spikemap(:)],ai,'linear'); % triplicate and interpolate
            ratemap = spikemap ./ dwellmap; % calculate HD firing rate

            % head direction analyses
            hd3n = ratemap ./ max(ratemap); % normalise cell hd
            hd3n = hd3n(:);            
            r = circ_r(ai,hd3n); % rayleigh vector length
            mx = rad2deg(ai(hd3n == max(hd3n))); % preferred angle (location of max frate)
            mn = rad2deg(circ_mean(ai,hd3n)); % mean angle
            sd = rad2deg(circ_std(ai,hd3n)); % std deviation angle
    end

    if isempty(r)
        r = NaN;
    end    
    if isempty(mx)
        mx = NaN;
    end
    if isempty(sd)
        sd = NaN;
    end
    if isempty(mn)
        mn = NaN;
    end




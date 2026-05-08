function [g,gdata] = get_grid_score(amap,binsize,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short desc.
% long description
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
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 21/08/17 created 
% version 1.1.0, Release 12/09/17 added Allen and Wills methods
% version 1.2.0, Release 13/09/17 added peak method
% version 1.2.1, Release 15/10/17 removed buggy peak method and corrected Allen method, added orientation
% version 2.0.0, Release 20/06/22 renamed get_grid_score, majorly overhauled 
% version 2.0.1, Release 20/06/22 added Krupic ellipticity measure
% version 2.0.2, Release 20/06/22 removed shared computations out of switch
% version 2.0.3, Release 20/06/22 simplified many of the methods
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
    p = inputParser;
    addRequired(p,'amap',@(x) isnumeric(x));  
        def_method = 'allen';    
        expectedmethods = {'allen','wills','langston','soman','mixed','brandon','sargolini','krupic','savelli'};  
    addRequired(p,'binsize',@(x) ~isempty(x) && ~all(isnan(x(:))));          
    addParameter(p,'method',def_method,@(x) any(validatestring(x,expectedmethods)));   
    addParameter(p,'figure',false,@(x) islogical(x));    
    parse(p,amap,binsize,varargin{:});
    config = p.Results;

    % preallocate grid values
    gdata.centroid = NaN(1,2); 
    gdata.diameter = NaN; 
    gdata.radius = NaN;     
    gdata.area = NaN;     
    gdata.convexarea = NaN; 
    gdata.field_orientation = NaN; 
    gdata.majaxislength = NaN;
    gdata.minaxislength = NaN; 
    gdata.elongation = NaN;   
    gdata.height = NaN; 
    gdata.width = NaN; 
    gdata.aspect_ratio = NaN;   
    gdata.grid_score = NaN;        
    gdata.grid_orientation = NaN;
    gdata.wavelength = NaN;      
    gdata.peaks_mask = NaN;
    g = NaN;
    if isempty(amap) || all(isnan(amap(:)))
        return
    end

    % based on the method chosen, retrieve the correct threshold for detecting the central peak
    peak_threshold = {'allen',0.1;'wills',0.3;'langston',0.3;'soman',NaN;'mixed',0.3;'sargolini',0.1;'brandon',0.1;'krupic',0;'savelli',0.1};
    threshold_now = peak_threshold{ find(ismember(peak_threshold(:,1),config.method),1) , 2};    
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FIELD STATS
%% Average firing field characteristics
% The central peak of the autocorrelogram can also provide information about the 'average'
% firing field size and shape. We have to find its position anyway, so we can just extract
% this information and add it to the output
    % find central peak
    im = imbinarize(amap,0.5); % threshold at 0.5
    bw = bwlabeln(im); % label peaks
    sts = regionprops('table',bw,'Centroid','EquivDiameter','Area','ConvexArea','Orientation','MajorAxisLength','MinorAxisLength','BoundingBox'); % get their properties
    sts = sts(sts.Area>=9,:); % remove small peaks
    sts.ds = pdist2(sts.Centroid,size(im,[2 1])./2); % find the distance of each peak from the center
    [~,i] = sortrows(sts.ds,'ascend'); % the central peak should be the first index in i (the 6 hexagonal peaks should also be indices 2:7)

    % accumulate central field data
    if isempty(sts) % if no peaks were found
        return
    end
    gdata.centroid = sts.Centroid(i(1),:); 
    gdata.diameter = sts.EquivDiameter(i(1),:) .* binsize; 
    gdata.radius = sts.EquivDiameter(i(1),:) ./ 2 .* binsize;           
    gdata.radius2 = sqrt( sts.Area(i) ) / pi;    
    gdata.area = sts.Area(i(1),:) .* (binsize^2);     
    gdata.convexarea = sts.ConvexArea(i(1),:) .* (binsize^2); 
    gdata.field_orientation = sts.Orientation(i(1),:); 
    gdata.majaxislength = sts.MajorAxisLength(i(1),:) .* binsize;
    gdata.minaxislength = sts.MinorAxisLength(i(1),:) .* binsize; 
    gdata.elongation = sts.MajorAxisLength(i(1),:) ./ sts.MinorAxisLength(i(1),:); 
    gdata.height = sts.BoundingBox(i(1),4) .* binsize;
    gdata.width = sts.BoundingBox(i(1),3) .* binsize;
    gdata.aspect_ratio = gdata.height / gdata.width;     

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% Calculate grid score
    % find central peaks
    im = imbinarize(amap,threshold_now); % threshold at threshold_now
    bw = bwlabeln(im); % label peaks
    sts = regionprops('table',bw,'Centroid','EquivDiameter','Area','ConvexArea','Orientation','MajorAxisLength','MinorAxisLength','BoundingBox'); % get their properties
    sts = sts(sts.Area>=9,:); % remove small peaks
    sts.ds = pdist2(sts.Centroid,size(im,[2 1])./2); % find the distance of each peak from the center
    [~,i] = sortrows(sts.ds,'ascend'); % the central peak should be the first index in i (the 6 hexagonal peaks should also be indices 2:7)

    %% Analyses common to all approaches
    npoints = sts.Centroid(i,:);
    dpoints = sts.ds(i,:);
    if strcmp(config.method,'brandon')
        npoints = npoints(dpoints < (dpoints(1)*1.5),:); % Brandon et al. specify that none of the six peaks could exceed 1.5 times the distance to the first         
    end    
    if size(sts.Centroid,1)<7 % if there are fewer peaks than we expect from a grid cell, just use whatever peaks there are
        hexagon_points = npoints(2:end,:);
        mds = mean(dpoints(2:end));                
    else
        hexagon_points = npoints(2:7,:);
        mds = mean(dpoints(2:7));                
    end

    hexagon_angles = atan2d(hexagon_points(:,2),hexagon_points(:,1)); % angle to each peak
          
    % calculate grid 'skew' which is estimated as the ellipticity of an ellipse fitted to the central 6 peaks
    % Method taken from:
    %   Krupic et al. (2015) Nature
    %   Grid cell symmetry is shaped by environmental geometry
    %   https://doi.org/10.1038/nature14153
    if size(hexagon_points,1)<5
        gdata.ellipticity = NaN;
        gdata.ellipse = [];
    else
        % a = fit_ellipse(hexagon_points(:,1),hexagon_points(:,2));
        % if ~isempty(a)
        %     gdata.ellipticity =  sqrt( 1 - (a.short_axis^2 / a.long_axis^2) );
        %     gdata.ellipse = a;
        % end
        gdata.ellipticity = NaN;
        gdata.ellipse = [];
    end
    
    % correct for Skew (still to add)
    % Savelli et al. (2017)    
    % If the elliptical index was >0.05, the rate map was ‘stretched’ along the direction of the 
    % shorter axis so as to correct the distortion. The autocorrelogram, the seven most central 
    % correlation fields, and their centers of mass were then recomputed from this rate map. 

    % distance transform of the autocorrelation
    dmat = zeros(size(amap)); % prepare a matrix of zeros, same size as autocorr
    dmat(round(size(amap,1)/2),round(size(amap,2)/2)) = 1; % centre blob equals 1       
    pixel_distances = bwdist(dmat); % pixel_distances contains distance of every pixel from the middle pixel
    
    %% Each approach in detail
    switch config.method       
        case {'wills','sargolini','krupic','allen','brandon','mixed','savelli'}
            %% Method used by:
            %   Sargolini et al. (2006)
            %   Conjunctive Representation of Position, Direction, and Velocity in Entorhinal Cortex
            %   https://doi.org/10.1126/science.1125572 \
            %   Note that Sargolini et al. define peaks as 100 or more contiguous pixels of 1.5 × 1.5 cm2 
            %   above a fixed threshold (in most cases r = 0.10) - the area cutoff seems far too high so I ignore it
            %   The correlation was confined to the area defined by a circle around the outermost peak of the six 
            %   peaks closest to the centre of the autocorrelation map
            %
            %   Wills et al. (2012) Front. Neural Circuits
            %   The abrupt development of adult-like grid cell firing in the medial entorhinal cortex
            %   https://doi.org/10.3389/fncir.2012.00021           
            %
            %   Krupic et al. (2012) Nature
            %   Grid cell symmetry is shaped by environmental geometry
            %   https://doi.org/10.1038/nature14153                  
            %   The six local peaks of the autocorrelogram were defined as the six local maxima with r > 0 closest 
            %   to the central peak (excluding the central peak itself). Gridness was calculated by defining a 
            %   mask on the spatial autocorrelogram centred on the central peak but excluding the peak itself 
            %   bounded by a circle defined by the mean distance from the centre to the closest peaks multiplied 
            %   by 2.5. This area was rotated in 30 degrees increments up to 150 degrees and for each rotation 
            %   the Pearson product-moment correlation coefficient was calculated against the un-rotated mask. 
            %   Gridness is then calculated taking the difference between the minimum correlation coefficient for 
            %   rotations of 60 and 120 degrees and the maximum correlation coefficient for rotations of 30, 90 
            %   and 150 degrees.     
            %
            %   Brandon et al. (2011)
            %   Reduction of Theta Rhythm Dissociates Grid Cell Spatial Periodicity from Directional Tuning
            %   https://doi.org/10.1126/science.1201652
            %
            %   Savelli et al. (2017)
            %   Framing of grid cells within and beyond navigation boundaries
            %   https://doi.org/10.7554/eLife.21354
            
            % calculate mean distance to closest blobs
            switch config.method
                case {'wills'}
                    ring_outer = mds*1.25;
                    inner_threshold = 0.5;
                    ring_inner = 0; % this will essentially be ignored                                          
                case {'sargolini'}
                    ring_outer = mds*1.25;
                    inner_threshold = 0.5;       
                    ring_inner = 0; % this will essentially be ignored                                                              
                case {'krupic'}
                    ring_outer = mds*2.5;
                    inner_threshold = 0.2; 
                    ring_inner = 0; % this will essentially be ignored                                                              
                case {'allen'}
                    ring_outer = mds*1.25;
                    inner_threshold = 0.2; % not described in paper              
                    ring_inner = mds*0.40;
                case {'brandon'}
                    ring_outer = mds+(sts.EquivDiameter(i(1),:)/2);
                    inner_threshold = 0.2; % not described in paper                                         
                    ring_inner = mds*0.50;  
                case {'savelli'}
                    ring_outer = mds+(sts.EquivDiameter(i(1),:)/2);
                    inner_threshold = 0.2; % not described in paper                                          
                    ring_inner = mds-(sts.EquivDiameter(i(1),:)/2);                    
                case {'mixed'}
                    min_rad = (sts.EquivDiameter(i(1),:)/2) * 1.5; % our minimum radii                    
                    ring_outer = mds+min_rad;
                    inner_threshold = 0.2; % not described in paper                                        
                    ring_inner = mds-min_rad;   
            end            

            % Cut to the central portion of the autocorrelation:   
            %   "Gridness was calculated by defining a mask on the spatial autocorrelogram centered on the 
            %   central peak, but excluding the peak itself (r > 0.5), bounded by a circle defined by the 
            %   mean distance from the center of the six closest peaks, multiplied by 1.25.
            ring_mask = pixel_distances<ring_outer & pixel_distances>ring_inner; % logical mask for central peak"
            imcent = amap;
            imcent(~ring_mask) = NaN; % only take data 
           
            % find central peak so we can remove it
            im = imbinarize(amap,inner_threshold); % threshold at inner_threshold
            bw = bwlabeln(im); % label peaks
            sts = regionprops('table',bw,'Centroid','Area','PixelIdxList'); % get their properties
            sts = sts(sts.Area>=9,:); % remove small peaks
            sts.ds = pdist2(sts.Centroid,size(im,[2 1])./2); % find the distance of each peak from the center
            [~,i] = sortrows(sts.ds,'ascend'); % the central peak should be the first index in i (the 6 hexagonal peaks should also be indices 2:7)
            imcent(sts.PixelIdxList{i(1)}) = NaN; % remove central peak

            % rotational correlation
            %   "Gridness was then expressed as the lowest correlation obtained for rotations of 60° and 120°, 
            %   versus the unrotated mask, minus the highest correlation obtained at 30°, 90°, or 150°."
            rs = NaN(5,1);
            as = 30:30:150;
            for a = 1:length(as)
                mrot = imrotate(imcent,as(a),'bilinear','crop');
                r = corrcoef(mrot(:),imcent(:),'rows','pairwise'); 
                rs(a) = r(1,2);
            end
 
            % save results
            gdata.grid_score = min([rs(2),rs(4)],[],'omitnan') - max([rs(1),rs(3),rs(5)],[],'omitnan'); 
            gdata.peaks_mask = ring_mask;
            gdata.grid_orientation = min(hexagon_angles,[],'omitnan');
            gdata.wavelength = mds;

        case {'langston'}
            %% Method based on that used by Langston et al. (2010) Development of the Spatial Representation System in the Rat
            % rings are cut from the autocorrelation at different distances from the centre and the standard rotation and correlation
            % is performed on each one. The highest of these grid scores is used, the radius of the ring is the grid spacing. A sine wave
            % is fitted to the values in the ring and this is used to estimate the positions of the fields. The grid orientation is the
            % angle from the centre to the first one of these fields, counter-clockwise. The grid field radius is taken as the radius of
            % the centre field.
            rad_steps = 1; % amount by which to increase radii each step
            min_radius = (sts.EquivDiameter(i(1),:)/2) * 1.5; % our minimum radii
            radius_width = min_radius;

            % calculate grid score on ever increasing radii
            ds = min_radius:rad_steps:((min(size(im))/2)-radius_width/2); % distances we will use
            if isempty(ds) % if the central peak is too large we can't form a ring
                return
            end

            gs = NaN(length(ds),1); % prepare an empty vector
            for dd = 1:length(ds) % for every distance we want to test
                ring_mask = pixel_distances>ds(dd)-radius_width & pixel_distances<ds(dd)+radius_width; % logical mask for central peak
                imcent = amap;
                imcent(~ring_mask) = NaN;                

                % rotational correlation
                rs = NaN(5,1); % prepare empty vector
                as = 30:30:150; % angles we want to test
                for a = 1:length(as) % for every angle
                    mrot = imrotate(imcent,as(a),'bilinear','crop'); % rotate the autocorr by this angle
                    r = corrcoef(mrot(:),imcent(:),'rows','pairwise'); % correlate the rotated map with the original
                    rs(a) = r(1,2); % collect the r value
                end
                gs(dd) = min([rs(2),rs(4)],[],'omitnan') - max([rs(1),rs(3),rs(5)],[],'omitnan'); % gs is the grid score    
            end
            [g,dindx] = max(gs,[],'omitnan'); % get the maximum grid score and find the distance associated with it

            % save results
            gdata.grid_score = g;        
            gdata.grid_orientation = min(hexagon_angles,[],'omitnan');
            gdata.wavelength = ds(dindx); % ds(dindx) is the distance with the best grid score         
            gdata.peaks_mask = pixel_distances>ds(dindx)-rad_width & pixel_distances<ds(dindx)+rad_width; % logical mask for central peak

        case {'koenig'}
            %% Method used by:
            %   Koenig et al. (2012) Science
            %   The Spatial Periodicity of Grid Cells Is Not Sustained During Reduced Theta Oscillations
            %   https://doi.org/10.1126/science.1201685         
            
            % An annulus that contains these surrounding peaks was identified by first
            % calculating the average autocorrelation for concentric circles from the center. The inner
            % boundary of the annulus was then set to the distance from the center at which the
            % concentric values were initially negative or at which they first reached a local minimum,
            % whichever was lower. This excluded the central peak. The outer boundary was set to the
            % second local minimum of the concentric values or at least 25 cm from the inner boundary.
            % The maps were then rotated in 30 degree steps, and the pixels that were contained in the
            % annulus were correlated. If 6-fold symmetry exists, the correlation values at 30, 90, and
            % 150 degrees are expected to be low and the values at 60 and 120 degrees are expected to
            % be high. ‘Gridness’ was defined as the average difference between the first and the
            % second set of values.      
            % As an alternate method, we performed the same calculation, but
            % averaged the correlation values of all pixels within the annulus that were along the same
            % angle from the center (i.e., along a ‘ray’) before performing the rotations. The correlation
            % between the pixel-wise method in (S11, S13) and the ray-wise method was 0.96, but the
            % latter method was more robust in identifying cells with small degrees of ellipticity along
            % one of the grid axis. The scores that were obtained with the latter method are thus
            % reported.          
            
        case {'soman'} 
            % Soman, K., Muralidharan, V., and Chakravarthy, V.S. (2017)
            % A Model of Multisensory Integration and its Influence on Hippocampal Spatial Cell Responses
            % IEEE Transactions on Cognitive and Developmental Systems

            % rotational correlation
            rs = NaN(5,1);
            as = 30:30:150;
            for a = 1:length(as)
                mrot = imrotate(im,as(a),'bilinear','crop');
                r = corrcoef(im(:),mrot(:),'rows','pairwise'); 
                rs(a) = r(1,2);
            end

            % save results
            gdata.grid_score = min([rs(2),rs(4)],[],'omitnan') - max([rs(1),rs(3),rs(5)],[],'omitnan');        
            gdata.grid_orientation = min(hexagon_angles,[],'omitnan');

    end
    g = gdata.grid_score;
    
    if isempty(gdata.grid_orientation)
        gdata.grid_orientation = NaN;
    end    
    if isempty(gdata.grid_score)
        gdata.grid_score = NaN;
    end  
    if isempty(gdata.wavelength)
        gdata.wavelength = NaN;
    end      
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUMMARY FIGURE
%% Produce a figure if requested
    % create figure if required
    if config.figure
        figure();
        imc = imagesc(amap);
        set(imc,'alphadata',gdata.peaks_mask.*0.8); hold on
        title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',gdata.grid_score,gdata.wavelength,gdata.radius,gdata.orientation));   
        xlabel('x (px)')
        ylabel('y (px)')
        caxis([0 max(imcent(:),[],'omitnan')])
        daspect([1 1 1]);
        axis off
    end     






























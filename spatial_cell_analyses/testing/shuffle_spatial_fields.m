function [shuffled_maps] = shuffle_spatial_fields(ratemap, n_shuffles)
% shuffle_spatial_fields Generates spatial null distributions via field deformation
%
% This function fragments a firing rate map into spatial segments, generates 
% target segments of identical sizes in a blank map, and bijectively maps the 
% original pixels to the new footprint using boundary control-point warping. 
% This preserves exact local firing rate statistics while disrupting global topography.
%
% USAGE
%
% [shuffled_maps] = shuffle_spatial_fields(ratemap, n_shuffles)
%
% INPUT
%
% 'ratemap'    - (double) A 2D matrix representing the smoothed spatial firing rate.
% 'n_shuffles' - (integer) The number of shuffled maps to generate.
%
% OUTPUT
%
% 'shuffled_maps' - (3D matrix) Dimensions [size(ratemap,1), size(ratemap,2), n_shuffles].
%                   Contains the generated null firing rate maps.
%
% NOTES
% 1. Target segments are grown using seeded dilation to guarantee exact pixel quotas.
% 2. Pixel assignment uses a greedy KNN approach on warped coordinates to avoid interpolation.
%
% EXAMPLE
% 
% % Generate 100 shuffled maps for a given place cell
% shuffled_maps = shuffle_spatial_fields(PC_ratemap, 100);
% 
% SEE ALSO get_synth_neurons

% HISTORY
%
% version 1.0.0, Release 04/09/26 Initial draft release
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        ratemap (:,:) double {mustBeNonempty}
        n_shuffles (1,1) double {mustBePositive, mustBeInteger} = 100
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    
    [rows, cols] = size(ratemap);
    shuffled_maps = zeros(rows, cols, n_shuffles);
    
    % Step 1: Detect original spatial segments (Module 1)
    orig_labels = segment_ratemap(ratemap);
    num_segments = max(orig_labels(:));
    
figure
tiledlayout('flow')
nexttile
imagesc(orig_labels)
% keyboard

    % Calculate exact pixel quotas required for each segment
    segment_quotas = zeros(1, num_segments);
    for s = 1:num_segments
        segment_quotas(s) = sum(orig_labels(:) == s);
    end

    for i = 1:n_shuffles
        % Step 2: Generate target footprint masks (Module 2)
        target_labels = generate_target_segments([rows, cols], segment_quotas);
        
        current_shuffled_map = zeros(rows, cols);
        


nexttile
imagesc(target_labels)
keyboard


        % Step 3: Warp and assign pixels per segment (Module 3)
        for s = 1:num_segments
            source_mask = (orig_labels == s);
            target_mask = (target_labels == s);
            
            mapped_values = warp_and_assign_pixels(ratemap, source_mask, target_mask);
            
            % Insert mapped values into the overall shuffled map
            current_shuffled_map(target_mask) = mapped_values;
        end
        
        shuffled_maps(:,:,i) = current_shuffled_map;
    end
end

%%%%%%%%%%%%%%%% HELPER FUNCTIONS
function [target_labels] = generate_target_segments(map_size, quotas)
% generate_target_segments Grows solid, contiguous spatial segments to exact quotas
%
% Generates target spatial footprints using synchronous, layer-by-layer distance 
% expansion from seeds. Guarantees that every region is a single solid, contiguous 
% blob while strictly enforcing exact pixel quotas.
%
% USAGE
%
% [target_labels] = generate_target_segments(map_size, quotas)
%
% INPUT
%
% 'map_size' - (double) 1x2 array specifying the [rows, cols] of the target map.
% 'quotas'   - (double) 1xN array specifying the exact number of pixels required 
%              for each segment.
%
% OUTPUT
%
% 'target_labels' - (double) A 2D matrix of size map_size containing integer 
%                   labels for each segment.
%
% NOTES
% 1. Guarantees zero interleaving; each region expands strictly outward from its seed.
%
% EXAMPLE
% 
% labels = generate_target_segments([100, 100], [500, 1000, 1500]);
% 
% SEE ALSO segment_ratemap, shuffle_spatial_fields

% HISTORY
%
% version 3.0.0, Release 04/09/26 Contiguous region growth release
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        map_size (1,2) double {mustBeNonnegative, mustBeInteger}
        quotas (1,:) double {mustBePositive, mustBeInteger}
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY

    r = map_size(1);
    c = map_size(2);
    num_seg = length(quotas);
    target_labels = zeros(r, c);
    
    total_pixels = r * c;
    if sum(quotas) > total_pixels
        error('Sum of segment quotas exceeds total available map pixels.');
    end

    % 1. Robust Seed Placement (Peak/Distance Sorting with Spacing)
    valid_mask = true(r, c); 
    dist_map = bwdist(~valid_mask) + rand(r, c) * 0.1;
    [~, sorted_idx] = sort(dist_map(:), 'descend');
    
    seeds_idx = zeros(1, num_seg);
    seeds_count = 0;
    min_separation = max(5, round(min(map_size) / (num_seg + 2)));
    
    for i = 1:length(sorted_idx)
        candidate = sorted_idx(i);
        [sy, sx] = ind2sub([r, c], candidate);
        
        if seeds_count == 0
            seeds_count = seeds_count + 1;
            seeds_idx(seeds_count) = candidate;
        else
            accepted_coords = [floor((seeds_idx(1:seeds_count)-1)/r)+1, mod(seeds_idx(1:seeds_count)-1, r)+1];
            dists_to_existing = sqrt((sx - accepted_coords(:,1)).^2 + (sy - accepted_coords(:,2)).^2);
            if all(dists_to_existing >= min_separation)
                seeds_count = seeds_count + 1;
                seeds_idx(seeds_count) = candidate;
            end
        end
        if seeds_count == num_seg; break; end
    end
    
    if seeds_count < num_seg
        remaining_needed = num_seg - seeds_count;
        unassigned_pool = setdiff(sorted_idx, seeds_idx(1:seeds_count));
        seeds_idx(seeds_count+1:num_seg) = unassigned_pool(1:remaining_needed);
    end

    for idx = 1:num_seg
        target_labels(seeds_idx(idx)) = idx;
    end

    % 2. Synchronous Front-Propagation (Guarantees Solid Contiguous Regions)
    % Expands all regions layer-by-layer simultaneously until quotas are hit.
    current_counts = ones(1, num_seg);
    se = strel('disk', 1); % Smooth circular structural element for organic growth
    
    while sum(current_counts < quotas) > 0
        active_segs = find(current_counts < quotas);
        active_segs = active_segs(randperm(length(active_segs))); % Randomize priority
        
        for i = 1:length(active_segs)
            seg = active_segs(i);
            if current_counts(seg) >= quotas(seg); continue; end
            
            mask = (target_labels == seg);
            dilated = imdilate(mask, se);
            new_pixels = dilated & (target_labels == 0);
            new_idx = find(new_pixels);
            
            if isempty(new_idx); continue; end
            
            space_needed = quotas(seg) - current_counts(seg);
            if length(new_idx) > space_needed
                % Pick a random subset of the frontier pixels to maintain smooth edges
                new_idx = new_idx(randperm(length(new_idx), space_needed));
            end
            
            target_labels(new_idx) = seg;
            current_counts(seg) = current_counts(seg) + length(new_idx);
        end
    end
    
    % 3. Final Safety Catch-All for Any Leftover Unassigned Pixels
    zero_idx = find(target_labels == 0);
    if ~isempty(zero_idx)
        for z = 1:length(zero_idx)
            [zy, zx] = ind2sub([r, c], zero_idx(z));
            neighbors = target_labels(max(1, zy-1):min(r, zy+1), max(1, zx-1):min(c, zx+1));
            valid_neighbors = neighbors(neighbors > 0);
            if ~isempty(valid_neighbors)
                target_labels(zero_idx(z)) = valid_neighbors(1);
            else
                target_labels(zero_idx(z)) = 1;
            end
        end
    end

end



    function [labels] = segment_ratemap(ratemap, opts)
    % segment_ratemap Segments a spatial firing rate map into cohesive fields
    %
    % Applies heavy Gaussian smoothing to the input ratemap to suppress noise and merge 
    % sub-peaks, followed by a watershed transform to identify basin boundaries. Ridge lines 
    % and segments smaller than a defined threshold are reassigned to their nearest valid 
    % neighbor to ensure complete and cohesive tiling of the spatial environment.
    %
    % USAGE
    %
    % [labels] = segment_ratemap(ratemap)
    % [labels] = segment_ratemap(ratemap, 'min_area', 10, 'smooth_sigma', 3)
    %
    % INPUT
    %
    % 'ratemap' - (double) A 2D matrix representing the spatial firing rate.
    % 'min_area' - (Name-Value) Scalar double specifying the minimum pixel count for a segment.
    %                           Default is 5.
    % 'smooth_sigma' - (Name-Value) Scalar double specifying the Gaussian smoothing standard deviation.
    %                               Default is 3.
    %
    % OUTPUT
    %
    % 'labels' - (double) A 2D matrix of identical size to the input, containing integer segment labels.
    %
    % NOTES
    % 1. The input ratemap is smoothed internally only for boundary detection. The output labels 
    %    are meant to mask the original un-smoothed (or lightly smoothed) ratemap values.
    % 2. Assumes a fully populated matrix (e.g., a square environment without NaNs).
    %
    % EXAMPLE
    % 
    % % Segment a firing rate map with a 10-pixel minimum area
    % [labels] = segment_ratemap(PC_ratemap, 'min_area', 10);
    % 
    % SEE ALSO shuffle_spatial_fields
    
    % HISTORY
    %
    % version 1.0.0, Release 04/09/26 Initial release
    %
    % AUTHOR 
    % Roddy M. Grieves
    % University of Glasgow, School of Psychology and Neuroscience
    % eMail: roddy.grieves@glasgow.ac.uk
    % Copyright 2026 Roddy M. Grieves
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
    %%%%%%%%%%%%%%%% ARGUMENT CHECK
        arguments
            ratemap (:,:) double {mustBeNonempty}
            opts.min_area (1,1) double {mustBePositive, mustBeInteger} = 5
            opts.smooth_sigma (1,1) double {mustBePositive} = 3
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    
        % 1. Heavy Smoothing (Isolates macro-structure)
        smoothed_map = imgaussfilt(ratemap, opts.smooth_sigma);
        
        % 2. Invert and Watershed
        inv_map = max(smoothed_map(:)) - smoothed_map;
        labels = watershed(inv_map);
        
        % 3. Minimum Area Filter
        num_segments = max(labels(:));
        for i = 1:num_segments
            if sum(labels(:) == i) < opts.min_area
                labels(labels == i) = 0; % Send small segments to 0 for reassignment
            end
        end
        
        % 4. Vectorized Nearest Neighbor Reassignment
        % Resolves both the natural ridge lines (0s) and the deleted micro-segments
        mask_valid = (labels > 0);
        
        if any(~mask_valid(:)) && any(mask_valid(:))
            [~, nearest_idx] = bwdist(mask_valid);
            labels(~mask_valid) = labels(nearest_idx(~mask_valid));
        end
        
        % 5. Re-index labels sequentially to ensure no missing integers
        unique_labels = unique(labels(:));
        for i = 1:length(unique_labels)
            labels(labels == unique_labels(i)) = i;
        end
    
    end


    function assigned_vals = warp_and_assign_pixels(full_map, s_mask, t_mask)
        % Warps source boundary to target boundary, then bijectively assigns pixels.
        
        % Get physical coordinates
        [sy, sx] = find(s_mask);
        [ty, tx] = find(t_mask);
        source_vals = full_map(s_mask);
        
        % Extract boundaries for control points
        s_bound = bwboundaries(s_mask);
        t_bound = bwboundaries(t_mask);
        
        % Fallback for extremely small/fragmented segments: bypass warp and just sort
        if isempty(s_bound) || length(s_bound{1}) < 10
            assigned_vals = sort(source_vals, 'descend'); % naive assignment
            return;
        end
        
        s_perim = s_bound{1};
        t_perim = t_bound{1};
        
        % Sample N control points evenly around the perimeter
        n_ctrl = min([10, size(s_perim,1), size(t_perim,1)]);
        s_idx = round(linspace(1, size(s_perim,1), n_ctrl));
        t_idx = round(linspace(1, size(t_perim,1), n_ctrl));
        
        s_ctrl = [s_perim(s_idx, 2), s_perim(s_idx, 1)]; % [x, y]
        t_ctrl = [t_perim(t_idx, 2), t_perim(t_idx, 1)];
        
        % Compute spatial transformation
        try
            tform = fitgeotrans(s_ctrl, t_ctrl, 'pwl');
            [warped_x, warped_y] = transformPointsForward(tform, sx, sy);
        catch
            % If geometry is collinear/degenerate, fall back to simple translation
            warped_x = sx - mean(sx) + mean(tx);
            warped_y = sy - mean(sy) + mean(ty);
        end
        
        % Vectorized Greedy KNN assignment to snap floating points to exact target grid
        ideal_coords = [warped_x, warped_y];
        target_coords = [tx, ty];
        
        assigned_vals = zeros(length(tx), 1);
        available_target_idx = (1:length(tx))';
        unassigned_source_idx = (1:length(sx))';
        
        while ~isempty(unassigned_source_idx)
            % Find nearest available grid bin for each unassigned warped point
            [idx_in_avail, dists] = knnsearch(target_coords(available_target_idx, :), ...
                                              ideal_coords(unassigned_source_idx, :));
                                          
            % Identify collisions (multiple source points claiming same target bin)
            unique_targets = unique(idx_in_avail);
            for k = 1:length(unique_targets)
                tgt = unique_targets(k);
                competitors = find(idx_in_avail == tgt);
                
                % Tie-break: keep the one with the smallest distance
                [~, winner_idx] = min(dists(competitors));
                true_winner_source = unassigned_source_idx(competitors(winner_idx));
                true_tgt_target = available_target_idx(tgt);
                
                % Assign value
                assigned_vals(tgt) = source_vals(true_winner_source); 
                
                % Remove from pools
                unassigned_source_idx(competitors(winner_idx)) = [];
                available_target_idx(tgt) = [];
                
                % Adjust indices for next loop iteration
                break; % Simplification for the draft: in a full implementation, we vectorize this removal
            end
        end
    end

function [grid_score, metrics, fft_mag] = compute_fourier_grid_score(rmap, rmset)
% compute_fourier_grid_score Computes Azimuthal Fourier grid score and spatial metrics
%
% Computes a 2D Fast Fourier Transform (FFT) magnitude from a firing rate map 
% using zero-padding, extracts grid spacing via radial profiling, computes an 
% azimuthal angular profile, and derives a Z-scored 6th-harmonic grid score 
% and global orientation analytically from Fourier phases following Krupic et al. (2012)[cite: 1].
%
% USAGE
%
% [grid_score, metrics, fft_mag] = compute_fourier_grid_score(rmap, rmset)
%
% INPUT
%
% 'rmap'  - (matrix) 2D firing rate map.
% 'rmset' - (struct) Rate map settings containing 'binsize' in mm.
%
% OUTPUT
%
% 'grid_score' - (scalar) Z-scored 6th-harmonic angular Fourier grid score.
% 'metrics'    - (struct) Extracted spatial metrics (spacing, wavelength, orientation, peaks).
% 'fft_mag'    - (matrix) 2D Fourier spectrogram magnitude (zero-padded/smoothed).
%
% NOTES
% 1. Uses azimuthal angular harmonics with true Z-score normalization and 
%    phase-based global orientation estimation, avoiding fragile local peak detection.
%
% EXAMPLE
% 
% [gs, metrics, fft_m] = compute_fourier_grid_score(rmap, rmset);
% 
% SEEALSO get_synth_neurons, rate_mapper

% HISTORY
%
% version 5.0.0, Release 24/08/26 Reverted to clean analysis pipeline with phase-based orientation
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        rmap (:,:) double
        rmset (1,1) struct
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % 1. Determine dimensions and bin size (mm to cm)
    [nrows, ncols] = size(rmap);
    binsize_cm = rmset.binsize / 10; % convert mm to cm

    % Fill NaN values with mean rate
    mean_rate = nanmean(rmap(:));
    rmap_filled = rmap;
    rmap_filled(isnan(rmap_filled)) = mean_rate;
    
    % 2. Zero-pad rate map for smooth Fourier spectrogram interpolation
    target_fft_size = 256;
    pad_r = max(0, target_fft_size - nrows);
    pad_c = max(0, target_fft_size - ncols);
    
    rmap_padded = padarray(rmap_filled, [floor(pad_r/2), floor(pad_c/2)], mean_rate, 'both');
    if size(rmap_padded, 1) < target_fft_size
        rmap_padded(end+1,:) = mean_rate;
    end
    if size(rmap_padded, 2) < target_fft_size
        rmap_padded(:,end+1) = mean_rate;
    end
    
    [p_rows, p_cols] = size(rmap_padded);
    padded_size_cm = p_cols * binsize_cm;
    
    % 3. 2D Fourier Transform & Shifted Magnitude
    rate_detrended = rmap_padded - mean(rmap_padded(:));
    F = fft2(rate_detrended);
    fft_mag = fftshift(abs(F));
    
    % 4. Dynamic Radial Profiling to find Grid Spacing / Wavelength
    cx = floor(p_cols / 2) + 1;
    cy = floor(p_rows / 2) + 1;
    
    [XGrid, YGrid] = meshgrid(1:p_cols, 1:p_rows);
    R_matrix = sqrt((XGrid - cx).^2 + (YGrid - cy).^2);
    
    max_possible_r = floor(min(cx, cy) * 0.9);
    min_r = 3; % Skip DC center
    r_bins = min_r:max_possible_r;
    radial_mean = zeros(size(r_bins));
    
    for idx = 1:length(r_bins)
        r_val = r_bins(idx);
        ring_mask = (R_matrix >= r_val - 0.5) & (R_matrix < r_val + 0.5);
        radial_mean(idx) = mean(fft_mag(ring_mask));
    end

    % Find peak radius corresponding to grid spacing
    [~, peak_idx] = max(radial_mean);
    r_peak = r_bins(peak_idx);
    
    wavelength = padded_size_cm / r_peak;
    spacing = wavelength;
    
    % 5. Compute 1D Angular Profile within Annular Band around r_peak
    ring_width = max(3, round(r_peak * 0.25));
    annulus_mask = (R_matrix >= r_peak - ring_width) & (R_matrix <= r_peak + ring_width);
    
    [Theta_matrix, ~] = cart2pol(XGrid - cx, YGrid - cy); % -pi to pi
    num_angle_bins = 360;
    angle_edges = linspace(-pi, pi, num_angle_bins + 1);
    angle_centers = angle_edges(1:end-1) + diff(angle_edges)/2;
    angular_profile = zeros(1, num_angle_bins);
    
    for a = 1:num_angle_bins
        bin_mask = annulus_mask & (Theta_matrix >= angle_edges(a)) & (Theta_matrix < angle_edges(a+1));
        angular_profile(a) = sum(fft_mag(bin_mask));
    end
    
    % 6. 1D FFT of Angular Profile for Z-Scored Grid Score & Phase Orientation
    ang_fft = fft(angular_profile - mean(angular_profile));
    ang_fft_mag = abs(ang_fft);
    ang_fft_phase = angle(ang_fft);
    
    % 6th harmonic is at index 7 (MATLAB 1-based indexing: index 1 is DC)
    h6_magnitude = ang_fft_mag(7);
    h6_phase = ang_fft_phase(7);
    
    % True Z-score calculation against background harmonics (excluding DC and 6th)
    background_harmonics = ang_fft_mag([2:6, 8:15]);
    grid_score = (h6_magnitude - mean(background_harmonics)) / (std(background_harmonics) + eps);
    
    % Analytical global orientation derived directly from 6th harmonic phase (modulo 60 degrees / pi/3)
    orientation = mod(h6_phase / 6, pi / 3);
    
    % 7. Populate metrics struct
    metrics.spacing = spacing;
    metrics.wavelength = wavelength;
    metrics.orientation = orientation;
    metrics.r_peak = r_peak;
    metrics.angular_profile = angular_profile;
    metrics.angle_centers = angle_centers;
    metrics.h6_magnitude = h6_magnitude;
    metrics.h6_phase = h6_phase;

%%%%%%%%%%%%%%%% HELPER FUNCTIONS
end
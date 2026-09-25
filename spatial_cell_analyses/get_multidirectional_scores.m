function [z_scores, peak_orientations] = get_multidirectional_scores(hd_rmap,target_k)
% get_multidirectional_scores Calculates Z-scored n-fold directional symmetry
%
% This function runs a 1D FFT over a circular head direction rate map to 
% quantify unidirectional (1-fold), bidirectional (2-fold), and 
% quad-directional (4-fold) tuning using a Z-score against background harmonics.
% It effectively quantifies multidirectional symmetry in spatial cells.
%
% USAGE
%
% [z_scores, peak_orientations] = get_multidirectional_scores(hd_rmap)
%
% INPUT
%
% 'hd_rmap'           - (numeric vector) 1D head direction firing rate map. 
%                       Assumes bins cover a full 360 degrees.
%
% OUTPUT
%
% 'z_scores'          - (numeric vector) [1 x 3] Z-scores for the 1st (uni), 
%                       2nd (bi), and 4th (quad) harmonics respectively.
% 'peak_orientations' - (numeric vector) [1 x 3] Orientation of the primary 
%                       peak (in radians) for the 1st, 2nd, and 4th harmonic 
%                       symmetries.
%
% NOTES
% 1. The Z-score is calculated by comparing the target harmonic magnitude against 
%    a background noise floor of other harmonics (excluding the target and DC).
%
% 2. The background is bounded to a maximum of the 15th harmonic. This prevents 
%    unphysiological, high-frequency bin-to-bin noise from artificially deflating 
%    the target harmonic's Z-score.
%
% EXAMPLE
% 
% % Generate a simulated noisy bidirectional tuning curve
% bins = linspace(0, 2*pi, 60);
% hd_rmap = cos(2 * bins) + randn(1, 60)*0.2 + 1.5;
% [z_scores, angles] = get_multidirectional_scores(hd_rmap);
% disp(z_scores); % The 2nd element (bidirectional) will be highly significant
% 
% SEE ALSO calc_spatial_metrics rate_mapper_hd

% HISTORY
%
% version 1.0.0, Release 26/08/26 Initial release
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        hd_rmap {mustBeNumeric, mustBeVector}
        target_k (1,:) double = [1,2,4]
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % Ensure it is a column vector and handle potential NaNs from unvisited bins
    hd_rmap = hd_rmap(:);
    hd_rmap(isnan(hd_rmap)) = 0; 
    
    % Number of bins and maximum harmonic to use for background (e.g., up to 15 cycles)
    N = length(hd_rmap);
    max_harmonic = min(15, floor(N/2) - 1); 
    
    % 1D FFT, centering the data first to suppress the DC component (index 1)
    hd_fft = fft(hd_rmap - mean(hd_rmap));
    hd_fft_mag = abs(hd_fft);
    hd_fft_phase = angle(hd_fft);
    
    % Define the target cycles we want to look for (1=Uni, 2=Bi, 4=Quad)
    % In MATLAB's 1-based indexing, harmonic k is at index k + 1
    % target_k = [1, 2, 4];
    
    % Preallocate outputs
    z_scores = NaN(1, 3);
    peak_orientations = NaN(1, 3);
    
    for i = 1:length(target_k)
        k = target_k(i);
        target_idx = k + 1; 
        
        % Extract magnitude and phase for the target harmonic
        target_mag = hd_fft_mag(target_idx);
        target_phase = hd_fft_phase(target_idx);
        
        % Define background harmonics (e.g., 1st through 15th harmonic)
        % We use index 2 to max_harmonic + 1 because index 1 is DC.
        bg_indices = 2:(max_harmonic + 1);
        
        % Remove the target harmonic from the background pool
        bg_indices(bg_indices == target_idx) = [];
        
        % Calculate Z-score against the background noise floor
        bg_mags = hd_fft_mag(bg_indices);
        z_scores(i) = (target_mag - mean(bg_mags)) / (std(bg_mags) + eps);
        
        % Calculate anatomical orientation of the symmetry (radians)
        peak_orientations(i) = mod(-target_phase / k, 2 * pi / k);
    end
end
function test_get_ahv_info()
% test_get_ahv_info Unit test and visualization for get_ahv_info
%
% Generates synthetic position timestamps and head direction angles using a 
% sinusoidal velocity profile, creates synthetic spikes, runs get_ahv_info, 
% and plots the resulting dwell, spike, and rate maps.
%
% USAGE
%
% test_get_ahv_info()
%
% INPUT
%
% None
%
% OUTPUT
%
% None
%
% NOTES
% 1. Uses a smooth sinusoidal head direction trajectory to test both 
%    positive and negative angular head velocities.
%
% EXAMPLE
% 
% test_get_ahv_info();
% 
% SEE ALSO get_ahv_info

% HISTORY
%
% version 1.0.0, Release 24/08/26 Initial release
%
% AUTHOR 
% Roddy M. Grieves
% University of Glasgow, School of Psychology and Neuroscience
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy M. Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        % No inputs required for this test script
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Generate Synthetic Data
    fprintf('Generating synthetic trajectory...\n');
    
    fs = 50;                  % Sampling rate (Hz)[cite: 1]
    dt = 1 / fs;
    duration_s = 30;          % 30 seconds of data
    pot = (0:dt:duration_s)'; % Position timestamps[cite: 1]
    
    % Create a head direction profile that oscillates smoothly (reversing direction)
    % Using a sine wave for position means velocity will be a cosine wave (smoothly varying AHV)
    oscillation_freq = 0.2;   % Hz
    max_hd_excursion = 4 * pi; % Total radians to sweep
    hd_rad = max_hd_excursion * sin(2 * pi * oscillation_freq * pot);
    poh = rad2deg(wrapToPi(hd_rad)); % Wrap between -180 and 180[cite: 1]

%%%%%%%%%%%%%%%% ARGUMENT CHECK (Synthetic Spikes)
    % Generate spikes preferentially when the animal is turning in one direction
    % Let's make the cell fire whenever the raw derivative of hd_rad is positive
    hd_velocity = gradient(hd_rad, dt);
    
    % Simple probabilistic or thresholded spiking: fire when velocity > threshold
    firing_prob = (hd_velocity > 1.5); 
    spindx = find(firing_prob);
    
    % Ensure spindx is non-empty
    if isempty(spindx)
        error('Synthetic spike generation failed: no spikes generated.');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY (Run Function)
    fprintf('Running get_ahv_info on synthetic data...\n');
    [scores, slopes, res] = get_ahv_info(pot, poh, spindx);

%%%%%%%%%%%%%%%% HELPER FUNCTIONS (Visualization)
    % Visualise the output maps
    figure('Name', 'AHV Unit Test Results', 'Color', 'w', 'Position', [100, 100, 900, 600]);
    
    subplot(2, 3, 1);
    plot(pot, poh, 'k');
    xlabel('Time (s)');
    ylabel('Head Direction (deg)');
    title('Synthetic Trajectory (POH)');
    
    subplot(2, 3, 2);
    plot(pot(spindx), poh(spindx), 'r.');
    xlabel('Time (s)');
    ylabel('Spike HD (deg)');
    title(sprintf('Spike Times (n = %d)', length(spindx)));
    
    subplot(2, 3, 3);
    bar(res.bin_centers, res.dwell_map, 'FaceColor', [0.5 0.5 0.5]);
    xlabel('AHV (deg/s)');
    ylabel('Dwell Time (s)');
    title('Dwell Map');
    
    subplot(2, 3, 4);
    bar(res.bin_centers, res.spike_map, 'FaceColor', [0 0 0]);
    xlabel('AHV (deg/s)');
    ylabel('Spike Count');
    title('Spike Map');
    
    subplot(2, 3, 5);
    bar(res.bin_centers, res.rate_map, 'FaceColor', [0 0.4 0.8]);
    xlabel('AHV (deg/s)');
    ylabel('Firing Rate (Hz)');
    title(sprintf('Rate Map (Left r: %.2f, Right r: %.2f)', scores(1), scores(2)));
    
    subplot(2, 3, 6);
    % Display summary metrics text box
    axis off;
    text(0.1, 0.7, sprintf('Scores [Left, Right]:\n  %.2f,  %.2f', scores(1), scores(2)), 'FontSize', 11, 'FontWeight', 'bold');
    text(0.1, 0.3, sprintf('Slopes [Left, Right]:\n  %.2f,  %.2f', slopes(1), slopes(2)), 'FontSize', 11);
    title('Tuning Metrics');

    fprintf('SUCCESS: Unit test executed and figures generated successfully.\n');
end
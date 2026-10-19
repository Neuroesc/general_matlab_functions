%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
% Define simulation parameters, epoch lengths, and map settings for the baseline test
n_hpc = 10;
n_grid = 0;
n_hd = 0;
n_noise = 10;

session_duration = 600; 
epoch_duration = 60; 

rmset.method = 'histogram';
rmset.binsize = 25; 
rmset.ssigma = 30; 
rmset.padding = 1;
rmset.mindwell = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% GENERATE SYNTHETIC DATA
% Simulate random walk with place cells and noise cells under baseline conditions (no noise/overdispersion)
[traj, cells] = get_synth_neurons(n_hpc, n_grid, n_hd, n_noise, ...
    'duration', session_duration, ...
    'noise_level', 0.5, ...
    'overdispersion', 0);

rmset.maplims = [min(traj.x(:),[],"all",'omitmissing') min(traj.y(:),[],"all",'omitmissing') max(traj.x(:),[],"all",'omitmissing') max(traj.y(:),[],"all",'omitmissing')].*10; 

%%%%%%%%%%%%%%%% DEFINE EPOCHS
% Split the continuous open-field session into discrete trials to satisfy ANOVA degrees of freedom
trial_starts = (0:epoch_duration:(session_duration - epoch_duration))';
epoch_trials = [trial_starts, trial_starts + epoch_duration];

%%%%%%%%%%%%%%%% PROCESS AND PLOT PLACE CELLS
% Extract spatial metrics for true place cells and plot their trajectories and spike locations
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_hpc
    [res, ~] = get_spatial_bootshuff(traj.x(:), traj.y(:), traj.t(:), rmset, ...
        'spt', cells.PC(i).t, ...
        'scores', {'spatial', 'anova'}, ...
        'trial', epoch_trials, 'eval_type',{'observed', 'shuffle', 'bootstrap'});

    % Extract SI metrics
    obs_si = res{'XY_si_spike', 'Observed'};
    p_obs_si = res{'XY_si_spike', 'p_value(right)_obs'};
    p_boot_si = res{'XY_si_spike', 'p_value(right)'};
    % keyboard

    % Extract ANOVA metrics
    obs_f = res{'ANOVA_F', 'Observed'};
    p_obs_f = res{'ANOVA_F', 'p_value(right)_obs'};
    p_boot_f = res{'ANOVA_F', 'p_value(right)'};

    % Plot trajectory and spikes
    nexttile;
    plot(traj.x, traj.y, 'Color', [0.8 0.8 0.8]); hold on;
    scatter(cells.PC(i).x, cells.PC(i).y, 10, 'r', 'filled');
    axis off; daspect([1 1 1]);

    % Format title with all 6 metrics
    title_str = sprintf('Obs SI: %.2f (p=%.3f) | Boot p=%.3f\nObs F: %.1f (p=%.3f) | Boot p=%.3f', ...
        obs_si, p_obs_si, p_boot_si, obs_f, p_obs_f, p_boot_f);
    title(title_str, 'FontSize', 8, 'FontWeight', 'normal');
end
sgtitle('Place Cells: Trajectories, Spikes, and Spatial Metrics', 'FontWeight', 'bold');

%%%%%%%%%%%%%%%% PROCESS AND PLOT NOISE CELLS
% Extract spatial metrics for non-spatial noise cells and plot their trajectories and spike locations
figure('Color', 'w', 'Position', [150, 150, 1200, 500]);
tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_noise
    [res, ~] = get_spatial_bootshuff(traj.x(:), traj.y(:), traj.t(:), rmset, ...
        'spt', cells.Noise(i).t, ...
        'scores', {'spatial', 'anova'}, ...
        'trial', epoch_trials);
    
    % Extract SI metrics
    obs_si = res{'XY_si_spike', 'Observed'};
    p_obs_si = res{'XY_si_spike', 'p_value(right)_obs'};
    p_boot_si = res{'XY_si_spike', 'p_value(right)'};
% keyboard    
    % Extract ANOVA metrics
    obs_f = res{'ANOVA_F', 'Observed'};
    p_obs_f = res{'ANOVA_F', 'p_value(right)_obs'};
    p_boot_f = res{'ANOVA_F', 'p_value(right)'};
    
    % Plot trajectory and spikes
    nexttile;
    plot(traj.x, traj.y, 'Color', [0.8 0.8 0.8]); hold on;
    scatter(cells.Noise(i).x, cells.Noise(i).y, 10, 'k', 'filled');
    axis off; daspect([1 1 1]);
    
    % Format title with all 6 metrics
    title_str = sprintf('Obs SI: %.2f (p=%.3f) | Boot p=%.3f\nObs F: %.1f (p=%.3f) | Boot p=%.3f', ...
        obs_si, p_obs_si, p_boot_si, obs_f, p_obs_f, p_boot_f);
    title(title_str, 'FontSize', 8, 'FontWeight', 'normal');
end
sgtitle('Noise Cells: Trajectories, Spikes, and Spatial Metrics', 'FontWeight', 'bold');
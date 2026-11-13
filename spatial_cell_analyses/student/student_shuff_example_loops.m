




n_place_cells = 2; % number of place cells we want per test, 100 would be good, less for testing
n_noise_cells = 2; % number of noise cells we want per test, 100 would be good, less for testing
session_duration = 1200; % length of simulated session in seconds, i.e. 20 mins = 1200 seconds
epoch_duration = 60; % length of "trials" in seconds, i.e. 1 min = 60 seconds

trial_starts = (0 : epoch_duration : (session_duration - epoch_duration))';
epoch_trials = [trial_starts, trial_starts + epoch_duration];

noise_levels = [0 2 4 6 8 10 20]; % values we want for noise, just an example for now
overdispersions = [0 0.2 0.4 0.6 0.8 1.0]; % values we want for dispersion, just an example for now

% DEFINE a structure containing the settings we want for generating firing rate maps
% Something like:
rmset = struct;
rmset.method = 'histogram'; 
rmset.binsize = 2.5; 
rmset.ssigma = 3.0; 
rmset.padding = 1;
rmset.mindwell = 0.02;

eval_types = ["observed", "shuffle", "bootstrap"];
metrics_to_run = {'spatial', 'anova'};
shuffle_iti = [2 2]; % shuffling iterations [bootstrap shuffle]

% The names and data types for our columns
varNames = {'noise', 'overdispersion', 'cell_type', ...
    'si_observed', 'si_bootstrapped', 'si_z_shuffle', 'si_z_bootshuff', ...
    'F_observed', 'F_bootstrapped', 'F_z_shuffle', 'F_z_bootshuff'};

% The total number of rows needed
n_cells = n_place_cells + n_noise_cells;
total_iterations = length(noise_levels) * length(overdispersions) * n_cells;

% Preallocate a matrix of NaNs and convert directly to a table
results = array2table(NaN(total_iterations,length(varNames),'single'), 'VariableNames', varNames);

row_counter = 1;
for n_idx = 1:length(noise_levels)
    for o_idx = 1:length(overdispersions)

        current_od = overdispersions(o_idx);
        current_noise = noise_levels(n_idx);
        fprintf('Processing OverD: %.2f, Noise: %.2f\n', current_od, current_noise);

        % CALL get_synth_neurons to generate Place Cells and Noise Cells 
        % USING the current_od and current_noise parameters
        [traj, cells] = get_synth_neurons(n_place_cells, 0, 0, n_noise_cells, 'duration', session_duration, ...
            'noise_level', current_noise, 'overdispersion', current_od); 
        rmset.maplims = [min(traj.epoly,[],1,'omitmissing') max(traj.epoly,[],1,'omitmissing')];

        % run through cells
        for cc = 1:2 % I'm going to do place then noise cells using the same loop, to save space
            if cc==1
                n_cells_now = n_place_cells;
            elseif cc==2
                n_cells_now = n_noise_cells;
            end

            for c_idx = 1:n_cells_now
                if cc==1
                    spt_now = cells.PC(c_idx).t;
                elseif cc==2
                    spt_now = cells.Noise(c_idx).t;
                end

                % CALL get_spatial_bootshuff on the place cells
                [res, ~] = get_spatial_bootshuff(traj.x(:), traj.y(:), traj.t(:), rmset, 'spt', spt_now, ...
                    'scores', metrics_to_run, 'trial', epoch_trials, 'eval_type', eval_types, 'iti', shuffle_iti);
    
                % enter data into our table
                res_si_row = find(ismember(res.Properties.RowNames,'XY_si_spike'));
                res_F_row = find(ismember(res.Properties.RowNames,'ANOVA_F'));
                res_Fp_row = find(ismember(res.Properties.RowNames,'ANOVA_p'));

                results.noise(row_counter,1) = current_noise;
                results.overdispersion(row_counter,1) = current_od;
                results.cell_type(row_counter,1) = cc; % 1 for place, 2 for noise cells

                results.si_observed(row_counter,1) = res.Observed(res_si_row); % observed 'actual' value
                results.si_bootstrapped(row_counter,1) = res.Boot_mean(res_si_row); % bootstrapped value (mean of bootstraps)
                results.si_z_shuffle(row_counter,1) = res.z_score_obs(res_si_row); % observed value z-scored relative to a shuffle
                results.si_z_bootshuff(row_counter,1) = res.z_score(res_si_row); % bootstrapped valye z-scored relative to a shuffle
                results.F_observed(row_counter,1) = res.Observed(res_F_row);
                results.Fp_observed(row_counter,1) = res.Observed(res_Fp_row);
                results.F_bootstrapped(row_counter,1) = res.Boot_mean(res_F_row);
                results.F_z_shuffle(row_counter,1) = res.z_score_obs(res_F_row);
                results.F_z_bootshuff(row_counter,1) = res.z_score(res_F_row);
                row_counter = row_counter + 1; % Increment our row counter with each cell
            end
        end
    end
end


% quick look at the output as an example
% isolate the true place cells
% A false negative can only occur if the cell is actually a place cell
pc_data = results(results.cell_type == 1, :);

% define the significance threshold
% Z = 1.645 roughly corresponds to p = 0.05 for a right-tailed test
z_thresh = 1.645; 

% identify false negatives (failed to reach significance)
is_false_neg = pc_data.si_z_shuffle < z_thresh;
% is_false_neg = pc_data.si_observed < 0.25;

% preallocate the heatmap matrix
fn_rate_matrix = NaN(length(overdispersions), length(noise_levels));

% 5. Bin the data and compute the average rate
for o_idx = 1:length(overdispersions)
    for n_idx = 1:length(noise_levels)

        % Extract current parameters
        curr_od = overdispersions(o_idx);
        curr_noise = noise_levels(n_idx);

        % Find all place cells belonging to this specific parameter bin
        bin_mask = (pc_data.overdispersion == curr_od) & (pc_data.noise == curr_noise);

        % Calculate the proportion of false negatives in this bin
        if any(bin_mask)
            fn_rate_matrix(o_idx, n_idx) = mean(is_false_neg(bin_mask));
        end
    end
end

% 6. Plot the resulting matrix as an image
figure('Color', 'w', 'Position', [100, 100, 500, 400]);
imagesc(noise_levels, overdispersions, fn_rate_matrix);
axis xy; % Flips the Y-axis so 0 is at the bottom
colorbar;
clim([0 1]); % Locks the colorbar limits from 0% to 100%

title('False Negative Rate (Observed SI vs Shuffle)');
xlabel('Noise Level');
ylabel('Overdispersion');










keyboard













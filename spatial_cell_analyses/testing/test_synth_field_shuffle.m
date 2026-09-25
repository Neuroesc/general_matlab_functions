function test_synth_grid_score()
% test_synth_grid_score Generates synthetic grid cells and evaluates Fourier grid score
%
% Runs get_synth_neurons to simulate animal trajectory and grid cell activity, 
% then applies compute_fourier_grid_score to evaluate the Fourier-based grid score 
% and extracted parameters against ground truth.
%
% USAGE
%
% test_synth_grid_score()
%
% INPUT
%
% None.
%
% OUTPUT
%
% None (generates diagnostic plots).
%
% NOTES
% 1. Uses get_synth_neurons to generate synthetic data.
%
% EXAMPLE
% 
% test_synth_grid_score();
% 
% SEE ALSO compute_fourier_grid_score, get_synth_neurons

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
    % No arguments required
    close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    fprintf('Generating synthetic neural data using get_synth_neurons[cite: 2]...\n');
    % [trajectory, cells] = get_synth_neurons(5, 5, 0, 'duration', 600, 'plot', false,'grid_spacing',[20 40]);
    % putvar(trajectory, cells);
    % boxsize = 100;
    % binsize = 2;
    trajectory = evalin('base','trajectory');
    cells = evalin('base','cells');

    % figure('Position', [100, 100, 1400, 700], 'Color', 'w');
    
    for ii = 1:2
        figure('Position', [100, 100, 1400, 700], 'Color', 'w');
    
        if ii==1
            num_grid_cells = length(cells.Grid);
        else
            num_grid_cells = length(cells.PC);
        end

        for i = 1:num_grid_cells
            if ii==1
                spikes_x = cells.Grid(i).x;
                spikes_y = cells.Grid(i).y;
            else
                spikes_x = cells.PC(i).x;
                spikes_y = cells.PC(i).y;
            end

            pos = [trajectory.x(:), trajectory.y(:)].*10;
            spk = [spikes_x(:), spikes_y(:)].*10;
    
            rmset = struct;
            rmset.method = 'histogram';
            rmset.binsize = 25;
            rmset.ssigma = 40;
            speedlift = [];
            [rmap,dmap,~,~,speedlift] = rate_mapper(pos,spk,rmset,speedlift);
            [shuffled_maps] = shuffle_spatial_fields(rmap, 10);
    
            subplot(4, num_grid_cells, i);
            plot(pos(:,1),pos(:,2),'k'); hold on;
            plot(spk(:,1),spk(:,2),'r.','MarkerSize',20);         
            xlabel('X (cm)'); ylabel('Y (cm)');
    
            subplot(4, num_grid_cells, i + num_grid_cells);
            imagesc( rmap);
            axis image xy
            colormap jet;
            % title(sprintf('Place Cell %d\nScore: %.3f | Spac: %.1fcm | Ori: %.1f', i, grid_score, metrics.spacing, metrics.orientation));
            xlabel('X (cm)'); ylabel('Y (cm)');
    
            % Plot shuffle
            subplot(4, num_grid_cells, i + num_grid_cells+1);
            imagesc( shuffled_maps(:,:,1));
            axis image xy
            colormap jet;
            % title(sprintf('Place Cell %d\nScore: %.3f | Spac: %.1fcm | Ori: %.1f', i, grid_score, metrics.spacing, metrics.orientation));
            xlabel('X (cm)'); ylabel('Y (cm)');

            % Plot shuffle
            subplot(4, num_grid_cells, i + num_grid_cells+2);
            imagesc( shuffled_maps(:,:,2));
            axis image xy
            colormap jet;
            % title(sprintf('Place Cell %d\nScore: %.3f | Spac: %.1fcm | Ori: %.1f', i, grid_score, metrics.spacing, metrics.orientation));
            xlabel('X (cm)'); ylabel('Y (cm)');


        end
        % sgtitle('Fourier-Based Grid Score Evaluation on Synthetic place cells');
    end



    fprintf('Evaluation complete. Displaying results.\n');

%%%%%%%%%%%%%%%% HELPER FUNCTIONS
end
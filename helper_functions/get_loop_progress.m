function loopout = get_loop_progress(loopin, res)
% get_loop_progress Fast console progress tracker for loops
%
% Displays a dynamically updating percentage of loop progress in the MATLAB 
% console. It calculates optimal breakpoint intervals to display progress 
% without significantly impacting loop execution time. 
%
% USAGE
%
% loopout = get_loop_progress(loopin) Initialize with default 1% resolution (100 updates)
%
% loopout = get_loop_progress(loopin, 10) Initialize with a custom 10% resolution (10 updates)
%
% INPUT
%
% 'loopin'     - (Numeric) On the first call, this is the total number of expected loops.
%                (Struct) On subsequent calls inside the loop, this is the output struct 
%                from the previous call.
%
% 'res'        - (Optional Numeric) The number of progress updates to display. 
%                Default is 100 (updates every 1%). Set to 10 for updates every 10%.
%
% OUTPUT
%
% 'loopout'    - Structure containing the loop state, counters, and timers. 
%                Must be passed back into the function on the next iteration.
%
% NOTES
% 1. To maximize speed, this function intentionally avoids MATLAB's modern 
%    argument validation blocks on the update cycle, relying on a simple 
%    integer comparison to bypass execution on 99% of calls.
% 2. Uses fprintf with backspaces (\b) to overwrite the console text in-place.
%
% EXAMPLE
% 
%   loops = 50000;
%   loopout = get_loop_progress(loops);
%   for ii = 1:loops
%       % ... core loop math here ...
%       loopout = get_loop_progress(loopout);    
%   end
%
% SEE ALSO waitbar, tic, toc

% HISTORY
%
% version 1.0.0, Release 27/06/19 Initial release
% version 2.0.0, Release 19/08/26 Performance overhaul, fprintf in-place updates
% version 2.0.1, Release 19/08/26 Fixed deletion bug
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK

    if isnumeric(loopin) && isscalar(loopin)
        % Initialisation (Runs once)
        if nargin < 2 || isempty(res)
            res = 100; % Default to 1% resolution
        end
        
        loopout = struct();
        loopout.numloops = loopin;
        loopout.loop = 1;
        loopout.start = tic;
        
        % Calculate exactly when to trigger a console update
        loopout.disploops = ceil(linspace(1, loopin, res + 1));
        loopout.disp_idx = 2; % Start looking for the second breakpoint
        loopout.next_disp = loopout.disploops(2);
        
        % Print the initial state
        fprintf('\tprocessing: %3d%%', 0);
        
    else
        % Update block (Runs every iteration - optimized for speed)
        loopout = loopin;
        
        % FAST BYPASS: If we haven't reached the next display threshold, 
        % increment and immediately return.
        if loopout.loop < loopout.next_disp
            loopout.loop = loopout.loop + 1;
            return;
        end
        
        % Display update
        idx = loopout.disp_idx;
        
        % Calculate current percentage mathematically
        percent = round(((idx - 1) / (length(loopout.disploops) - 1)) * 100);
        
        % Delete the previous 4 characters (e.g., "  0%") and print the new one
        fprintf(repmat('\b', 1, 4));
        fprintf('%3d%%', percent);
        
        % Update state for the next milestone
        loopout.disp_idx = idx + 1;
        
        if loopout.disp_idx <= length(loopout.disploops)
            loopout.next_disp = loopout.disploops(loopout.disp_idx);
        else
            loopout.next_disp = inf; % Prevent further triggers
        end
        
        % Final completion stats on the very last loop
        if loopout.loop == loopout.numloops
            elapsed = toc(loopout.start);
            fprintf(' (%.2fs, %.4fs/loop)\n', elapsed, elapsed / loopout.numloops);
        end
        
        loopout.loop = loopout.loop + 1;  
    end
end
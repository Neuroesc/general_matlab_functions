function check_toolboxes(filename)
    % 1. Get required products for the specified file
    try
        [~, pList] = matlab.codetools.requiredFilesAndProducts(filename);
        required_toolboxes = {pList.Name}';
    catch
        error('Could not analyze %s. Make sure it is in your path.', filename);
    end
    
    % 2. Get currently installed products
    installed = ver;
    installed_toolboxes = {installed.Name}';
    
    % 3. Find the difference
    missing_toolboxes = setdiff(required_toolboxes, installed_toolboxes);
    
    % 4. Display results
    fprintf('\n--- Toolbox Check for %s ---\n', filename);
    if isempty(missing_toolboxes)
        fprintf('All required toolboxes are currently installed.\n\n');
    else
        fprintf('WARNING: You are missing the following toolboxes:\n');
        for i = 1:length(missing_toolboxes)
            fprintf('  - %s\n', missing_toolboxes{i});
        end
        fprintf('\n');
    end
end
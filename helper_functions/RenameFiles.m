function RenameFiles(old_name,new_name)                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Renaming files...');
disp(sprintf('\t...old filename is: %s',old_name));
disp(sprintf('\t...new filename is: %s',new_name));

files_to_change = dir;
files_to_change_mat = struct2cell(files_to_change);
files_to_change_names = files_to_change_mat(1,:);

disp(sprintf('\t...potential files found: %s',num2str(length(files_to_change_names)-2)));

num = 0;                                                                                                                                % start a counter
for name = 1:length(files_to_change_names)
	curr_file = files_to_change_names{name};
	dots = strfind(curr_file,'.');

	if ~isempty(dots)
		first_dot = dots(1);
		curr_file_start = curr_file(1:first_dot-1);

		if strcmp(curr_file_start,old_name);
			try
				file1 = curr_file;                                                                                               % create the file's current name							
				file2 = [curr_file, '.temp'];	                                                                                                % create a temporary name					
				file3 = strrep(curr_file,old_name,new_name);                                                                                        % create the new lower case name							
				movefile(file1,file2);	                                                                                                % move from current name to the temporary one						
				movefile(file2,file3);	                                                                                                % move from the temporary name to the new lower case one						
				num = num + 1;	 
			catch
				continue
			end
		end % if curr_file_start == old_name
	end % if length(dots)>0
end % for name = 1:length(files_to_change_names)

if (num == 0)                                                                                                                           % if no filenames were changed
        disp(sprintf('\t...no files needed to be changed'));                                                                  % make a message
else                                                                                                                                    % if some names were changed
	disp(sprintf('\t...%d files were changed', num));                                                                     % make a message
end % if (num == 0)

disp(sprintf('\t...done'));                                                                                                             % display 'done' - indented in the command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

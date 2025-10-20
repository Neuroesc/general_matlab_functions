function package_function(fname,dout)
% package_function pack a function for easy sharing
% copy a function and all of its dependencies into a single directory
%
% USAGE
%
% package_function(fname,dout)
%
% INPUT
%
% 'fname' - String, name of the function to package
%
% 'dout' - String, optional, directory in which to create the package, default
%           is current directory

% HISTORY
%
% version 1.0.0, Release 18/08/16 created 
% version 1.0.0, Release 18/08/16 added compression of resulting directory
% version 1.0.0, Release 02/09/16 improved comments
% version 1.0.0, Release 08/10/24 simplified
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2024 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    % if no output directory was given use the default
    if ~exist('dout','var') || isempty(dout)
        dout = [pwd '\' fname '_packed']; 
    end 
    [~,~,~] = mkdir(dout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % list dependencies of parent function
    [flist,~] = matlab.codetools.requiredFilesAndProducts(fname);

    % run through dependencies and copy them to the directory
    for ff = 1:length(flist)
        fnow = flist{ff};
        [~,nme,~] = fileparts(fnow);
        if strcmp(nme,fname)
            copyfile(fnow,dout,'f');        
        else
            [~,~,~] = mkdir([dout '\dependencies']);
            copyfile(fnow,[dout '\dependencies'],'f');
        end
    end 

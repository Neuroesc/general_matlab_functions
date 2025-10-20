function out = m_template(in,varargin)
% name short description
% longer description
%
% USAGE
%
% out = name(in) process with default settings
%
% out = name(in,optional) process using optional argument 1
%
% out = rate_mapper(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'in' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% OUTPUT
%
% 'out' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% NOTES
% 1. 
%
% 2. 
%
% EXAMPLE
% 
% SEE ALSO NAME, NAME

% HISTORY
%
% version 1.0.0, Release 00/00/25 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    % Parse inputs
    p = inputParser;
    addRequired(p,'in',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addOptional(p,'param2',def_param1,@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'param3',def_param2,@(x) isnumeric(x) && isscalar(x));   
    parse(p,in,varargin{:});
    config = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY

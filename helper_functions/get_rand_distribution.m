function r = get_rand_distribution(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%FUNCTION  short desc.
%
% USAGE:
%     out = template(in) process with default settings
% 
%     out = template(in,optional1) process using optional argument 1
% 
%     out = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%     of the process
% 
%     Parameters include:
% 
%     'param1'          -   (default = X) Scalar value, parameter to do something
% 
%     'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%     in    - input as a vector
%
% OUTPUT:
%     out   - output as a vector
%
% EXAMPLES:
%
%     % run function using default values
%     out = template(in,varargin)
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2020 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings
    def_dist                = 'normal';
    possible_dists          = {'normal','gaussian','poisson','beta','uniform','flat','gamma'}; 
    def_num                 = [100 1];        
    def_max                 = inf;        
    def_min                 = -inf;   
    
    % distribution specific values
    def_alpha               = 1;        
    def_beta                = 1;        
    def_lambda              = 1;        
    def_mean                = 1;        
    def_std                 = 1;        
       
%% Parse inputs
    p = inputParser;
    addParameter(p,'dist',def_dist,@(x) any(validatestring(x,possible_dists)));    
    addParameter(p,'num',def_num,@(x) isnumeric(x));   
    addParameter(p,'max',def_max,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'min',def_min,@(x) isnumeric(x) && isscalar(x));   
    
    % distribution specific values    
    addParameter(p,'alpha',def_alpha,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'beta',def_beta,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'lambda',def_lambda,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'mean',def_mean,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'std',def_std,@(x) isnumeric(x) && isscalar(x));   
    parse(p,varargin{:});

%% Retrieve parameters 
    config = p.Results;
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
% initial values
dnum = config.num;
if numel(dnum)==1
    dnum = [dnum 1];
end
dmin = config.min;
dmax = config.max;
        
%% situation where we want single output
if dmin == dmax % if the minimum and maximum are the same, just use the value with no distribution
    r = ones(dnum) .* dmin;
    return
end

%% generate the dstributions
switch config.dist
    case {'uniform','flat'}
        r = dmin + (dmax-dmin).*rand(dnum); % Uniform

    case {'normal','gaussian'}
        dmu = config.mean;
        dstd = config.std;

        count = 0;
        r = NaN(dnum);
        while count ~= numel(r(:))
            val = normrnd(dmu,dstd,1,1); % Gaussian
            if val >= dmin && val <= dmax
                count = count+1;
                r(count) = val;
            end 
        end 

    case {'beta'}
        dalpha = config.alpha;
        dbeta = config.beta;

        count = 0;
        r = NaN(dnum);
        while count ~= numel(r(:))
            val = betarnd(dalpha,dbeta,1,1) .* dmax; % Beta
            if val >= dmin && val <= dmax
                count = count+1;
                r(count) = val;
            end 
        end 

    case {'poisson'}
        dlambda = config.lambda;

        count = 0;
        r = NaN(dnum);
        while count ~= numel(r(:))
            val = poissrnd(dlambda,1,1); % Poisson
            if val >= dmin && val <= dmax
                count = count+1;
                r(count) = val;
            end 
        end 

    case {'gamma'}
        dalpha = config.alpha;
        dbeta = config.beta;

        count = 0;
        r = NaN(dnum);
        while count ~= numel(r(:))
            val = gamrnd(dalpha,dbeta,1,1); % gamma
            if val >= dmin && val <= dmax
                count = count+1;
                r(count) = val;
            end 
        end 

end



































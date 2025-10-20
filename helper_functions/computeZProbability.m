function [z,p,q] = computeZProbability(obs,shuff)
% computeZProbability probability (p-value) from shuffle
% given observed valuea and corresponding shuffle distributions this function
% calculates the z-scored values and two-tailed probability (p-value) or, if
% required the left and right one-tailed probability
%
% USAGE
%
% [z,p,q] = computeZProbability(obs,shuff) 
%
% INPUT
%
% 'obs' - [1xM], actual, observed values
%
% 'shuff' - [NxM] shuffle results, columns correspond to the values in obs,
%           rows correspond to shuffles. Empty rows or missing values should be 
%           NaNs.
%
% OUTPUT
%
% 'z' - [1xM], observed values, z-scored relative to corresponding shuffles
%
% 'p' - [1xM], p-values, likelihood of observing obs value(s) given the shuffle(s)
%       two-tailed
%
% 'q' - [2xM], p-values, likelihood of observing obs value(s) given the shuffle(s)
%       one-tailed, first row are left tailed probabilities, bottom row are
%       right-tailed probabilities.
%
% NOTES
% 1. This function calculates the corresponding p-value using the normal CPDF 
% (cumulative probability density function of the normal distribution). If you 
% want to make directional inferences (say something about the direction or 
% sign of the effect), use the one-tailed p-value, which corresponds to a one-sided 
% composite null hypothesis. If the direction of the effect does not matter, 
% use the two-tailed p-value, which corresponds to a point null hypothesis.
%
% 2. The p-value is used in the context of a Null-Hypothesis statistical test 
% (NHST) and it is the probability of observing the result which was observed, 
% or a more extreme one, assuming the null hypothesis is true 1. It is 
% essentially the proportional area of the Z distribution cut off at the point 
% of the observed Z score.
%
% EXAMPLE
% 
% [z,p,q] = computeZProbability(2,normrnd(0,1,1000,1));
% The observed value provided here is 2
% The shuffle provided is a normal distribution with a mean of 0 and a SD of 1
% Thus, the observed value should be around the 95th percentile of the shuffle
% The z result should always be roughly 1.96, which is the correct distance from
% the mean in units of standard deviation
% The p result should always be roughly 0.05, which is the correct probability
% given this z-value.
% Changing the observed value to something larger will decrease the p-value,
% smaller will increase the p-value.
% two-tailed p-values should be around 1-0.05/2 and 0.05/2 respectively
% 
% SEE ALSO normcdf, mean, std

% HISTORY
%
% version 1.0.0, Release 08/10/24 Initial release, previous copy lost
% version 1.0.1, Release 05/04/25 Renamed from stats_z_probability, improved comments
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    if size(obs,2) ~= size(shuff,2)
        error('inputs do not match in size: %d observed values but %d shuffles given',size(obs,2),size(shuff,2))
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % z-score values relative to shuffles
    z = (obs - mean(shuff,1,'omitmissing')) ./ std(shuff,[],1,'omitmissing'); % z-scored observed values

    % calculate p-values based on these
    p = 2*normcdf(-abs(z)); % two-tailed probability of z-scores
    p(p==0) = 1 ./ numel(shuff);    
    q = [normcdf(z); normcdf(z,'upper')]; % left and right one-tailed probability of z-scores
    q(q==0) = 1 ./ numel(shuff);    




























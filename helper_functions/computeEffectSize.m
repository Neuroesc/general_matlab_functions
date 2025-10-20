function e = computeEffectSize(v1,v2,varargin)
%computeEffectSize calculate effect sizes between two groups
% Calculate effect sizes like Cohen's d and Hedge's g between
% two groups.
%
% USAGE
%
%   e = computeEffectSize(v1,v2); calculates effect sizes for
%   the difference between group v1 and v2 which are output in
%   structure e.
%
% INPUT
%
%   'v1' - Mx1, numeric, group 1 data
%
%   'v2' - Mx1, numeric, group 1 data
%
% OUTPUT
%
%   outputs included in e output structure:
%
%   'mean_1'        -   Float, mean of group 1
%
%   'mean_2'        -   Float, mean of group 2
%
%   'sd_1'          -   Float, sample standard deviation of group 1
%
%   'sd_2'          -   Float, sample standard deviation of group 2
%
%   'num_1'         -   Scalar, number of values in group 1
%
%   'num_2'         -   Scalar, number of values in group 2
%
%   'mean_diff'     -   Float, difference in means between the groups
%
%   'df'            -   Scalar, degrees of freedom or the number of 
%                       values in group 1 and 2 minus 2 (Bessel correction)
%
%   'sd_pooled'     -   Float, pooled standard deviation of both groups
%                       with Bessel correction
%
%   'cohen_d'       -   Cohen's d, Cohen (1988), see note 1
%
%   'hedge_g'       -   Hedge's g, or the corrected Cohen's d for small
%                       group sizes, Hedges and Olkin (1985), see note 1
%
%   'glass_delta'   -   Glass's delta, Lakens (2013), see note 2
%
%   'cliff_delta'   -   Cliff's delta, Hess and Kromrey (2004)
%
%   'probability_of_superiority'   -   Probability of superiority
%                                      Grissom and Kim (2012)
%
%
%   Class Support
%   -------------
%   The input vectors v1 and v2 in must be a real, non-sparse matrix of
%   the following classes: uint8, int8, uint16, int16, uint32, int32,
%   single or double.
%
%   Notes
%   -----
%   1. There is a lot of confusion surrounding what should and should not be
%       called Cohen's d or Hedge's g, for simplicity I have followed the narrative
%       in Lakens (2013) Calculating and reporting effect sizes to facilitate 
%       cumulative science: a practical primer for t-tests and ANOVAs
%       https://doi.org/10.3389%2Ffpsyg.2013.00863
%       Where Cohen's d includes the Bessel correction and Hedge's g refers to
%       the unbiased correction for Cohen's d.
%
%   2. In Glass's delta Cohen's d is calculated usind the standard deviation of
%       the control group as a reference. In this function the second input, v2
%       is taken as the control group.
%
%   Example
%   ---------
%
%   v1 = normrnd(1,1,100,1); % input vector 1
%   v2 = normrnd(10,1,100,1); % input vector 2
%   e = computeEffectSize(v1,v2); % effect sizes
% 
%   See also permutationZTest, computeZProbability

% HISTORY:
% version 1.0.0, Release 18/11/23 Initial release
% version 1.0.1, Release 05/04/25 Renamed from stats_effect_size, improved comments
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    p = inputParser;
    addRequired(p,'v1',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addRequired(p,'v2',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    % addParameter(p,'test','',@(x) ischar(x) || isstring(x));   
    parse(p,v1,v2,varargin{:});
    config = p.Results;
    v1 = config.v1(:);
    v2 = config.v2(:);

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% General values
    e.mean_1 = mean(v1,'all','omitnan'); % mean of group 1
    e.mean_2 = mean(v2,'all','omitnan'); % mean of group 2
    e.sd_1 = std(v1,[],'all','omitnan'); % sample SD of group 1
    e.sd_2 = std(v2,[],'all','omitnan'); % sample SD of group 2
    e.num_1 = numel(v1); % N of group 1
    e.num_2 = numel(v2); % N of group 2
    e.mean_diff = e.mean_1 - e.mean_2; % difference in group means

    % Lakens (2013) Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs
    % https://doi.org/10.3389%2Ffpsyg.2013.00863
    % equation 1 (denominator)    
    e.df = e.num_1+e.num_2-2; % group degrees of freedom (with Bessel's correction for bias in the estimation of the population variance)
    e.sd_pooled = sqrt( (((e.num_1-1)*e.sd_1^2) + ((e.num_2-1)*e.sd_2^2)) / e.df ); % pooled standard deviation

%% Cohen's d
    % Lakens (2013) Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs
    % https://doi.org/10.3389%2Ffpsyg.2013.00863
    % equation 1
    % Cohen's d is used to describe the standardized mean difference of an effect. 
    % This value can be used to compare effects across studies, even when the dependent 
    % variables are measured in different ways, for example when one study uses 
    % 7-point scales to measure dependent variables, while the other study uses 
    % 9-point scales, or even when completely different measures are used, such as 
    % when one study uses self-report measures, and another study used physiological 
    % measurements. It ranges from 0 to infinity. Cohen (1988) refers to the standardized 
    % mean difference between two groups of independent observations for the sample as 
    % d which is given by:
    e.cohen_d = (e.mean_1 - e.mean_2) / e.sd_pooled; % Cohen's d
    % In this formula, the numerator is the difference between means of the two groups of 
    % observations. The denominator is the pooled standard deviation. Remember that the 
    % standard deviation is calculated from the differences between each individual 
    % observation and the mean for the group. These differences are squared to prevent 
    % the positive and negative values from cancelling each other out, and summed (also 
    % referred to as the sum of squares). This value is divided by the number of observations 
    % minus one, which is Bessel's correction for bias in the estimation of the population
    % variance, and finally the square root is taken. This correction for bias in the
    % sample estimate of the population variance is based on the least squares estimator
    % (see McGrath and Meyer, 2006).

%% Hedge's g  
    % Lakens (2013) Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs
    % https://doi.org/10.3389%2Ffpsyg.2013.00863
    % equation 4
    % As mentioned earlier, the formula for Cohen's d, which is based on sample 
    % averages gives a biased estimate of the population effect size (Hedges and 
    % Olkin, 1985), especially for small samples (n < 20). Therefore, Cohen's d 
    % is sometimes referred to as the uncorrected effect size. The corrected effect 
    % size, or Hedges's g (which is unbiased, see Cumming, 2012), is:
    e.hedge_g = e.cohen_d * (1 - (3/(4*(e.num_1+e.num_2)-9))); % Hedge's g
    % Although the difference between Hedges's g and Cohen's d is very 
    % small, especially in sample sizes above 20 (Kline, 2004), it is preferable (and 
    % just as easy) to report Hedges's g. There are also bootstrapping procedures to 
    % calculate Cohen's d when the data are not normally distributed, which can provide 
    % a less biased point estimate (Kelley, 2005). 

%% Glass's delta
    % Lakens (2013) Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs
    % https://doi.org/10.3389%2Ffpsyg.2013.00863
    % Table 1
    % Glass's delta is defined as the mean difference between the experimental and control group 
    % divided by the standard deviation of the control (second) group.
    e.glass_delta = (e.mean_1 - e.mean_2) / e.sd_2; % Glass's delta

%% Cliff's delta
    % Cliff, N. (1996). Ordinal methods for behavioral data analysis
    % Hess and Kromrey (2004) Robust Confidence Intervals for Effect Sizes: A Comparative Study 
    % of Cohen s d and Cliff’s Delta Under Non-normality and Heterogeneous Variances
    % equation 3
    % The sample estimate of this statistic, Cliff’s ˆδ , is obtained by comparing each of the scores 
    % in one group to each of the scores in the other. The calculation of this sample statistic is 
    % given by: 
    e.cliff_delta = ( sum(v1>v2','all') - sum(v1<v2','all') ) / (e.num_1*e.num_2); 
    % The non-parametric nature of Cliff s δ reduces the influence of such characteristics as distribution
    % shape, differences in dispersion and extreme values. The statistic relies on what Cliff refers
    % to as a dominance analysis, a concept referring to the degree to which one sample overlaps another: 
    % the greater the overlap (i.e., the lower the dominance), the less difference between the groups. 
    % Unlike Cohen s d, Cliff’s effect size is bounded. An effect size of 1.0 or -1.0 indicates the 
    % absence of overlap between the two groups whereas a 0.0 indicates no overlap and the group 
    % distributions are equivalent.

%% Probability of superiority
    % Ruscio and Gera (2013) Generalizations and Extensions of the Probability of Superiority Effect Size Estimator
    % https://doi.org/10.3389%2Ffpsyg.2013.00863
    % equation 1
    e.probability_of_superiority = ( sum(v1>v2','all') + (sum(v1==v2','all')*0.5) ) / (e.num_1*e.num_2); 
    % In other words, one makes all pairwise
    % comparisons between members of Group 1 and members of Group 2, tallies the
    % number of times that the former scores higher than the latter (or credits this
    % as .5 if they are tied), and divides by the total number of comparisons. Values
    % can range from .00 (all members of Group 2 score higher than all members
    % of Group 1) through .50 (equal probability that members of either group score
    % higher, more formally known as stochastic equality) to 1.00 (all members of
    % Group 1 score higher than all members of Group 2). A serves as an estimator
    % of  D Pr.Y1 > Y2/, which Grissom and Kim (2012) called the “probability of
    % superiority” (p. 149), and its calculation is simple and intuitive.












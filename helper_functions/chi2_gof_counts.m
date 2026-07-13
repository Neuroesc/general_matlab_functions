function [p, x, e] = chi2_gof_counts(N1, N2, options)
% chi2_gof_counts Chi-Square goodness-of-fit test for two categories
%
% Performs a two-sided Chi-Square Goodness-of-Fit test on summarized count data 
% (e.g., N of neuron 1 and N of neuron 2) to determine if the observed proportions 
% significantly differ from an expected null distribution (default 50/50).
%
% USAGE
%
% [p, chi2] = chi2_gof_counts(N1, N2) process with default 50/50 expected probability.
%
% [p, chi2] = chi2_gof_counts(N1, N2, 'p1', val) process with a custom Name-Value pair for the null hypothesis.
%
% INPUT
%
% 'N1' - Scalar, positive integer that specifies the observed count of category 1.
%       units are raw counts (i.e. neurons).
%       No default value; must be specified.
%
% 'N2' - Scalar, positive integer that specifies the observed count of category 2.
%       units are raw counts (i.e. neurons).
%       No default value; must be specified.
%
% 'p1' - (Optional Name-Value) Scalar, numeric value between 0 and 1.
%       Specifies the expected probability (chance rate) of an observation 
%       belonging to category 1 under the null hypothesis.
%       Default value is 0.5.
%
% OUTPUT
%
% 'p' - Scalar numeric value between 0 and 1.
%       Represents the probability of observing this large of a difference 
%       purely by chance (the p-value). Two-sided.
%
% 'x' - Scalar numeric value.
%       The calculated Chi-Square test statistic.
%
% 'e' - 1x2, numeric, expected values
%
% NOTES
% 1. This function uses the direct mathematical calculation for the Chi-Square 
%    statistic and chi2cdf to derive the p-value. This avoids binning issues 
%    associated with passing pre-summarized data to MATLAB's built-in chi2gof.
%
% 2. The expected probability for category B is automatically calculated as 
%    (1 - p1).
%
% EXAMPLE
% 
% % Test if 40 type 1 and 60 type 2 neurons deviate from a 50/50 split:
% [p, x] = chi2_gof_counts(40, 60);
%
% % Test if 40 type 1 and 60 type 2 deviate from an expected 25/75 split:
% [p, x] = chi2_gof_counts(40, 60, 'p1', 0.25);
%
% SEE ALSO chi2gof, chi2cdf, binomtest
%
% HISTORY
%
% version 1.0.0, Release 13/07/26 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    % Validate inputs
    arguments
        N1 (1,1) {mustBeInteger, mustBeNonnegative}
        N2 (1,1) {mustBeInteger, mustBeNonnegative}
        options.p1 (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(options.p1, 0), mustBeLessThanOrEqual(options.p1, 1)} = 0.5
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % setup observed counts and totals
    observed = [N1, N2];
    total_count = N1 + N2;
    
    % handle edge case where total counts are zero
    if total_count == 0
        p = NaN;
        x = NaN;
        warning('Total count is zero. Cannot compute chi-square statistic.');
        return;
    end

    % calculate expected counts based on the null hypothesis
    expectedA = total_count * options.p1;
    expectedB = total_count * (1 - options.p1);
    expected = [expectedA, expectedB];

    % handle edge case where an expected count is zero
    if any(expected == 0)
        p = NaN;
        x = NaN;
        warning('One of the expected counts is zero, check p1');
        return
    end

    % calculate Chi-Square test statistic: sum((O - E)^2 / E)
    x = sum((observed - expected).^2 ./ expected);
    e = expected;

    % calculate p-value (chi-square with 1 d.f.)
    p = 1 - chi2cdf(x, 1);
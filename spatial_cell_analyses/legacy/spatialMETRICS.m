 function smetric = spatialMETRICS(rmap,dmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spatialMETRICS  calculates spatial measures
% Takes a ratemap and a dwellmap and calculates a number of spatial measures such as spatial information content and sparsity
% Some of these are adapted to work for n dimensional data
%
% USAGE:
%         smetric = spatialMETRICS(rmap,dmap)
%
% INPUT:
%         rmap - firing rate map (Hz)
%         dmap - dwell time map (s)
%
% OUTPUT:
%    smetric - structure containing spatial measures
%        smetric.shannons_entropy = Shannon's entropy, defined as -sum(p.*log2(p)) (equiavlent to Matlab function entropy), lower entropy means the ratemap is more complex and thus has a higher information content: https://en.wikipedia.org/wiki/Entropy_(information_theory)
%        smetric.cross_entropy = cross entropy, defined as -sum(p.*log2(q)), higher cross-entropy means the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Cross_entropy
%        smetric.mutual_info = mutual information, how much does knowing about the ratemap tell you about the dwellmap
%        smetric.kl_divergence = Kullback–Leibler divergence, defined as sum(p.*log2(p)-log2(q)), higher KLdivergence means the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
%        smetric.kl_divergence_symmetric = symmetrised Kullback–Leibler divergence, higher values indicate the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence#Symmetrised_divergence
%        smetric.jensen_shannon_divergence = derived from KLd, bounded (0,1] the square root of this is the Jensen-Shannon distance metric, higher values indicate the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
%        smetric.jensen_shannon_distance = the square root of js divergence: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
%        smetric.spatial_information = Skagg's spatial information content (bits per second), higher values indicate the ratemap is less related to the dwellmap: http://markus.lab.uconn.edu/wp-content/uploads/sites/1005/2015/01/Markus1994-Spatial-Info.pdf
%        smetric.spatial_information_perspike = Skaggs spatial information divided by the overall firing rate (bits per spike)
%        smetric.sparsity = the relative proportion of the maze on which the cell fired, i.e.  a sparsity score of 0.10 indicates that the cell fired on 10% of the maze surface: http://markus.lab.uconn.edu/wp-content/uploads/sites/1005/2015/01/Markus1994-Spatial-Info.pdf
%        smetric.mean_method_focus = mean method focus measure (Helmli, 2001), gives a measure of the local smoothness of 9 pixel neighborhoods, high values mean the ratemap is rougher (this would correspond with being more focused, which is what the measure was for originally), I think this might substitute spatial coherence
%        smetric.signal_to_noise = the max ratemap value / the mean ratemap value, called selectivity by Skaggs et al. (1996)
%
% EXAMPLES:
%
% See also: ENTROPY MUTUALINFORMATION IMHIST

% HISTORY:
% version 1.0.0, Release 2017/04/04: Initial release created to contain all spatial measures in one place
% version 1.0.1, Release 2017/04/04: vectorised and checked spatial information and sparsity
% version 1.0.2, Release 2017/05/04: added mean method focus (spatial coherance)
% version 2.0.0, Release 2018/06/15: renamed spatialMETRICS, updated spatial information and sparsity, added detailed comments
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL ARGUMENTS CHECK
% convert to greyscale
rmap8 = im2uint8(rmap);
dmap8 = im2uint8(dmap);

% calculate histogram counts
p = imhist(rmap8(:));
q = imhist(dmap8(:));

% normalize distributions so that their sum is one.
p = p ./ numel(rmap8);
q = q ./ numel(dmap8);
p = p + eps; % add smallest possible amount to avoid log of 0
q = q + eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
%% Mutual information
%     Intuitively, mutual information measures the information that X and Y share: It measures how much knowing one of these variables reduces 
%     uncertainty about the other. For example, if X and Y are independent, then knowing X does not give any information about Y and vice versa, 
%     so their mutual information is zero. At the other extreme, if X is a deterministic function of Y and Y is a deterministic function of 
%     X then all information conveyed by X is shared with Y: knowing X determines the value of Y and vice versa. As a result, in this case the 
%     mutual information is the same as the uncertainty contained in Y (or X) alone, namely the entropy of Y (or X). Moreover, this mutual 
%     information is the same as the entropy of X and as the entropy of Y. (A very special case of this is when X and Y are the same random variable.)
%     https://en.wikipedia.org/wiki/Mutual_information
mi = MutualInformation(rmap8(:),dmap8(:));

%% Shannon's entropy (rmap only)
%     help for matlab function entropy which does the same thing:
%     E = entropy(I) returns E, a scalar value representing the entropy of an
%     intensity image.  Entropy is a statistical measure of randomness that can be
%     used to characterize the texture of the input image.  Entropy is defined as
%     -sum(p.*log2(p)) where p contains the histogram counts returned from IMHIST.
se = -sum(p.*log2(p));

%% Cross entropy
%     In information theory, the cross entropy between two probability distributions {\displaystyle p} p and {\displaystyle q} q over the same underlying 
%     set of events measures the average number of bits needed to identify an event drawn from the set, if a coding scheme is used that is optimized for 
%     an "unnatural" probability distribution {\displaystyle q} q, rather than the "true" distribution {\displaystyle p} p.
ce0 = p .* log2(q); % using log base 2 means the output is measured in bits
ce = -nansum(ce0(:));

%% Kullback–Leibler divergence
%     In mathematical statistics, the Kullback–Leibler divergence (also called relative entropy) is a measure of how one probability distribution diverges from 
%     a second, expected probability distribution. In the simple case, a Kullback–Leibler divergence of 0 indicates that we can expect similar, if not the same, 
%     behavior of two different distributions, while a Kullback–Leibler divergence of 1 indicates that the two distributions behave in such a different 
%     manner that the expectation given the first distribution approaches zero. In simplified terms, it is a measure of surprise.
%     https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
kld = sum(p .* (log2(p) - log2(q)));

%% Symmetric Kullback–Leibler divergence
%     This quantity has sometimes been used for feature selection in classification problems, where P and Q are the conditional pdfs of a feature under two different classes.
KL1 = sum(p .* (log2(p) - log2(q)));
KL2 = sum(q .* (log2(q) - log2(p)));
kls = (KL1 + KL2) / 2;

%% Jensen-Shannon divergence
%     In probability theory and statistics, the Jensen–Shannon divergence is a method of measuring the similarity between two probability distributions.
%     It is based on the Kullback–Leibler divergence, with some notable (and useful) differences, including that it is symmetric and it is always a finite value. 
%     The square root of the Jensen–Shannon divergence is a metric often referred to as Jensen-Shannon distance.
%     https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
logQvect = log2((q + p) / 2);
jsd = 0.5 * (sum(p .* (log2(p) - logQvect)) + sum(q .* (log2(q) - logQvect)));
jsdist = sqrt(jsd);

%% Skaggs spatial information content (bits per second)
%     Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
%     https://onlinelibrary.wiley.com/doi/epdf/10.1002/%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K
%     Markus et al. (1994) Spatial information content and reliability of hippocampal CA1 neurons: effects of visual input
%     https://onlinelibrary.wiley.com/doi/epdf/10.1002/hipo.450040404
pi = dmap ./ nansum(dmap(:)); % dwell time probability
ro = nansum(rmap(:) .* pi(:)); % overall firing rate
si = nansum(pi(:) .* (rmap(:)./ro) .* log2(rmap(:)./ro)); 

%% Skaggs spatial information content (bits per spike)
%     From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
%     The information rate given by formula (1) is measured in bits per second. If it is
%     divided by the overall mean ring rate of the cell (expressed in spikes per second),
%     then a different kind of information rate is obtained, in units of bits per spike|let us
%     call it the information per spike. This is a measure of the specificity of the cell: the
%     more grandmotherish" the cell, the more information per spike. 
sis = si ./ ro;

%% Sparsity
%     From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
%     The sparsity measure is  an adaptation  to space of  a formula invented by Treves and Rolls (1 99 1); the adaptation measures the 
%     fraction of the environment  in which a cell  is active. Intuitively, a sparsity of, say, 0.1 means that the place field of the cell 
%     occupies 1/10 of the area the rat traverses.
sp = (nansum(pi(:).*rmap(:)).^2) ./ (nansum(pi(:).*(rmap(:).^2)));

%% Mean method focus measure (Helmli, 2001)
%     This measure is being used in place of spatial coherence as it works in n-dims and is basically the same idea
%     From Helmli and Scherer (2001) Adaptive Shape from Focus with an Error Estimation in Light Microscopy 
%     https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=938626
%     If the image is becoming sharper the variance Of the grey Of that Scene is getting higher. The ratio Of the mean grey value to the center 
%     grey value in the neighborhood can be taken as a new focus measure. The ratio results in one if there is a constant grey value or there 
%     is no texture present. If high variations are present the ratio is different from one.
meanf = ones([5 5])./25;
U = imfilter(rmap8,meanf,'replicate');
R1 = U ./ rmap8;
R1(rmap8 == 0) = 1;
index = (U > rmap8);
fm = 1./R1;
fm(index) = R1(index);
fm = nanmean(fm(:));

%% Spatial Coherence (Cacucci et al. 2007)
%     https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007
%     The spatial coherence for each firing rate map was computed as the mean correlation between the firing rate of each bin with the 
%     aggregate rate of the 24 nearest bins.
meanf = ones([5 5]);
meanf(3,3) = 0;
U = imfilter(rmap,meanf,'replicate');
scohe = corr(rmap(:),U(:),'rows','pairwise','type','Pearson');

%% Signal to noise ratio / selectivity
%     Called selectivity by Skaggs et al. (1996)
%     The selectivity measure is  equal to the spatial maximum  fir- ing rate  divided  by  the mean  firing  rate  of the cell. The more tightly 
%     concentrated the cell’s activity, the higher  the selectivity. A cell  with  no spatial tuning at all will  have  a selectivity of  1; 
%     there is  in principle  no upper  limit. A similar measure was  used by Barnes et al. (1983), except that the “out-of-field” firing rate was 
%     used instead of the mean rate. The present definition is prefer- able because it does  not depend  on identifying  a  “place field,” and because 
%     it is much less sensitive to noise. 
snr = nanmax(rmap(:)) ./ nanmean(rmap(:));

%% Accumulate data in structure
smetric = struct; % structure to hold data
smetric.shannons_entropy = se;
smetric.cross_entropy = ce;
smetric.mutual_info = mi;
smetric.kl_divergence = kld;
smetric.kl_divergence_symmetric = kls;
smetric.jensen_shannon_divergence = jsd;
smetric.jensen_shannon_distance = jsdist;
smetric.spatial_information = si;
smetric.spatial_information_perspike = sis;
smetric.sparsity = sp;
smetric.mean_method_focus = fm;
smetric.signal_to_noise = snr;
smetric.spatial_coherence = scohe;









































function spat = getSPATinfo(rmap,dmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes a firing rate map and a dwell time map and calculates a number of spatial measures which might be useful
%   For many of these the rmap and dmap are greyscaled and transformed into probability distributions (p and q respectively) 
%   using imhist
%   The most popular of the included measures are probably spatial information and sparsity
%   These measures should all work equally well for 3D data
%   spat = getSPATinfo(rmap,dmap)
%
%%%%%%%% Inputs
%   rmap = firing rate map
%   dmap = dwell time map (same size as ratemap)
%
%%%%%%%% Outputs
%   spat = structure containing spatial measures:
%        spat.shannons_entropy = Shannon's entropy, defined as -sum(p.*log2(p)) (equiavlent to Matlab function entropy), lower entropy means the ratemap is more complex and thus has a higher information content: https://en.wikipedia.org/wiki/Entropy_(information_theory)
%        spat.cross_entropy = cross entropy, defined as -sum(p.*log2(q)), higher cross-entropy means the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Cross_entropy
%        spat.kl_divergence = Kullback–Leibler divergence, defined as sum(p.*log2(p)-log2(q)), higher KLdivergence means the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
%        spat.kl_divergence_symmetric = symmetrised Kullback–Leibler divergence, higher values indicate the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence#Symmetrised_divergence
%        spat.jensen_shannon_divergence = derived from KLd, bounded (0,1] the square root of this is the Jensen-Shannon distance metric, higher values indicate the ratemap is less related to the dwellmap: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
%        spat.spatial_information = Skagg's spatial information content (bits), higher values indicate the ratemap is less related to the dwellmap: http://markus.lab.uconn.edu/wp-content/uploads/sites/1005/2015/01/Markus1994-Spatial-Info.pdf
%        spat.sparsity = the relative proportion of the maze on which the cell fired, i.e.  a sparsity score of 0.10 indicates that the cell fired on 10% of the maze surface: http://markus.lab.uconn.edu/wp-content/uploads/sites/1005/2015/01/Markus1994-Spatial-Info.pdf
%        spat.mean_method_focus = mean method focus measure (Helmli, 2001), gives a measure of the local smoothness of 9 pixel neighborhoods, high values mean the ratemap is rougher (this would correspond with being more focused, which is what the measure was for originally), I think this might substitute spatial coherence
%        spat.signal_to_noise = the max ratemap value / the mean ratemap value
%
%%%%%%%% Comments
%   04/04/17 created to contain all spatial measures in one place
%   04/04/17 vectorised and checked spatial information and sparsity
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% convert to greyscale
rmap8 = im2uint8(rmap);
dmap8 = im2uint8(dmap);

% calculate histogram counts
p = imhist(rmap8(:));
q = imhist(dmap8(:));

% normalize distributions so that their sum is one.
p = p ./ numel(rmap8);
q = q ./ numel(dmap8);
p = p + eps; % add smallest possible amount to avoid log of 0
q = q + eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate spatial mesures, vectorally, on any dimension maps
% calculate mutual information
%     Intuitively, mutual information measures the information that X and Y share: It measures how much knowing one of these variables reduces 
%     uncertainty about the other. For example, if X and Y are independent, then knowing X does not give any information about Y and vice versa, 
%     so their mutual information is zero. At the other extreme, if X is a deterministic function of Y and Y is a deterministic function of 
%     X then all information conveyed by X is shared with Y: knowing X determines the value of Y and vice versa. As a result, in this case the 
%     mutual information is the same as the uncertainty contained in Y (or X) alone, namely the entropy of Y (or X). Moreover, this mutual 
%     information is the same as the entropy of X and as the entropy of Y. (A very special case of this is when X and Y are the same random variable.)
%     https://en.wikipedia.org/wiki/Mutual_information
mi = MutualInformation(rmap8(:),dmap8(:));

%% perform shuffle to determine probability of MI and obtain a scaled (z-scored) value
%     To correct the estimates of spatial information, we first shuffled the position of firing rates on each trial. This approach avoids any firing 
%     preference across trials and in the mean rate. We then computed the Ispike, Isec and MI metrics using the shuffled rates. This nprocedure was 
%     repeated 100 times and used to build a surrogate distribution for each metric. The information values were then corrected by expressing them as 
%     z-scores of the surrogate distribution. 
%     Souza et al. (2017) On information metrics for spatial coding http://dx.doi.org/10.1101/189084
%     Basically, these guys found that mutual information is better reflected with bayesian decoding than traditional spatial information content
%     Thus is should be better than SI. This works for 2D and 3D. The MI score is also scaled using a shuffled distribution.

% calculate Shannon's entropy (rmap only)
se0 = p .* log2(p); % using log base 2 means the output is measured in bits
se0(p == 0) = [];
se = -nansum(se0(:)); % could also just use Matlab function 'entropy'

% calculate cross entropy
ce0 = p .* log2(q); % using log base 2 means the output is measured in bits
ce = -nansum(ce0(:));

% calculate the Kullback–Leibler divergence
kld = sum(p .* (log2(p) - log2(q)));

% calculate the symmetric Kullback–Leibler divergence
KL1 = sum(p .* (log2(p) - log2(q)));
KL2 = sum(q .* (log2(q) - log2(p)));
kls = (KL1 + KL2) / 2;

% calculate the Jensen-Shannon divergence
logQvect = log2((q + p) / 2);
jsd = 0.5 * (sum(p .* (log2(p) - logQvect)) + sum(q .* (log2(q) - logQvect)));

% calculate Skaggs spatial information content (bits per second)
pi = dmap ./ nansum(dmap(:)); % dwell time probability
ro = nansum(rmap(:) .* pi(:)); % overall firing rate
si = nansum(pi(:) .* (rmap(:)./ro) .* log2(rmap(:)./ro)); % Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences

% calculate Skaggs spatial information content (bits per spike)
sis = si ./ ro;

% calculate sparsity
sp = (nansum(pi(:).*rmap(:)).^2) ./ (nansum(pi(:).*(rmap(:).^2)));

% calculate mean method focus measure (Helmli, 2001)
meanf = ones([3 3 3])./27;
U = imfilter(rmap8,meanf,'replicate');
R1 = U ./ rmap8;
R1(rmap8 == 0) = 1;
index = (U > rmap8);
fm = 1./R1;
fm(index) = R1(index);
fm = nanmean(fm(:));

% calculate signal to noise ratio / selectivity
snr = nanmax(rmap(:)) ./ nanmean(rmap(:));

%% accumulate data in structure
spat = struct; % structure to hold data
spat.shannons_entropy = se;
spat.cross_entropy = ce;
spat.mutual_info = mi;
spat.kl_divergence = kld;
spat.kl_divergence_symmetric = kls;
spat.jensen_shannon_divergence = jsd;
spat.spatial_information = si;
spat.spatial_information_perspike = sis;
spat.sparsity = sp;
spat.mean_method_focus = fm;
spat.signal_to_noise = snr;
























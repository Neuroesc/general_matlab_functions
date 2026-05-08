function m = get_spatial_info(dmap,rmap,opts)
% get_spatial_info calculate various spatial information measures
% Function calculates the most commonly used spatial information theoretic
% measures such as Skaggs spatial information content and Mutual Information,
% can operate on 1D (e.g. Head direction), 2D or 3D maps. The function only
% needs a firing rate map (in Hz) and a dwell time map (in seconds).
%
% USAGE
%
% m = get_spatial_info(dmap,rmap) process with default settings
%
% m = get_spatial_info(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'dmap' - MxN numeric matrix, required, dwell time map
%       units are in seconds.
%
% 'rmap' - MxN numeric matrix, required, firing rate map
%       units are in Hz (spikes per second).
%
% 'metrics' - Cell array of strings, optional
%       Choices include: 
%           {'all'} - compute all spatial metrics
%           {'spatial_info'} - compute Skaggs spatial information content, both bits per second and bits per spike
%           {'sparsity'} - Skaggs sparsity measure
%           {'mutual_info'} - Mutual Information
%           {'entropy'} - Shannon's entropy and Cross entropy
%           {'kld'} - Kullback–Leibler divergence, symmetric Kullback–Leibler divergence and Jensen-Shannon divergence
%           {'spatial_coherence'} - Cacucci et al.'s (2007) spatial coherence measure
%           {'signal_to_noise'} - Signal to noise ratio (max / mean)
%       default is {'all'}
%
% OUTPUT
%
% 'm' - Structure, fields will depend on the outputs requested:
%       m.skaggs_si_bits_per_sec = Skaggs (1996) spatial information content (bits per second)
%       m.skaggs_si_bits_per_spike = Skaggs (1996) spatial information content (bits per spike)
%       m.sparsity = Skaggs (1996) sparsity
%       m.mutual_info = Mutual information
%       m.shannon_entropy = Shannon's entropy
%       m.cross_entropy = Cross entropy
%       m.kldivergence = Kullback–Leibler divergence
%       m.kldivergence_symmetric Symmetric Kullback–Leibler divergence
%       m.jsdivergence = Jensen-Shannon divergence
%       m.spatial_coherence = Cacucci et al. (2007) spatial coherence
%       m.snr = signal to noise (max / mean)
%
% NOTES
% 1. See the section of code for each metric for more detail about where the
% method was adapted from and how it works
%
% 2. Some metrics do not work on the original dwell and rate maps and instead
% use a greyscale version of binned values, you can see this operation at the
% top of the code
%
% EXAMPLE
%
% % synthetic ratemap and dwellmap
% mapSize = 100;
% numFields = 2;
% sigma = 8;
% rmap = zeros(mapSize);
% rmap( randperm(mapSize * mapSize, numFields) ) = 1;
% rmap = imgaussfilt(rmap,sigma,'FilterSize',2*ceil(4*sigma)+1,'Padding',0);
% rmap = (rmap / max(rmap(:))) * 15;
% dmap = rand(mapSize);
% dmap = imgaussfilt(dmap, 4);
% dmap = (dmap / sum(dmap(:))) * 600;
% 
% % get the spatial metrics
% m = get_spatial_info(dmap,rmap,'metrics',{'all'});
% 
% % plot results
% figure; tiledlayout; nexttile;
% imagesc(rmap); 
% title(sprintf('Ratemap\nSpatial info: %.2fb/s',m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
% axis image; colorbar;
% 
% % synthetic ratemap and dwellmap
% numFields = 8;
% sigma = 16;
% rmap = zeros(mapSize);
% rmap( randperm(mapSize * mapSize, numFields) ) = 1;
% rmap = imgaussfilt(rmap,sigma,'FilterSize',2*ceil(4*sigma)+1,'Padding',0);
% rmap = (rmap / max(rmap(:))) * 15;
% 
% % get the spatial metrics
% m = get_spatial_info(dmap,rmap,'metrics',{'spatial_info','sparsity'});
% 
% % plot results
% nexttile;
% imagesc(rmap); 
% title(sprintf('Ratemap\nSpatial info: %.2fb/s',m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
% axis image; colorbar;
% 
% nexttile
% imagesc(dmap); 
% axis image; colorbar;
% title(sprintf('Dwellmap'),'FontWeight','normal'); 
% 
% SEE ALSO 

% HISTORY
%
% version 1.0.0, Release 08/05/26 Initial release, adapted from spatialMETRICS
% version 1.0.1, Release 08/05/26 Cleaned up comments, inputs, code layout
% version 1.1.0, Release 08/05/26 Greatly improved choice of spatial metrics
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        dmap
        rmap
        opts.metrics (1,:) string {mustBeMember(opts.metrics, ["all", "spatial_info", "sparsity", "mutual_info", "entropy", "kld", "spatial_coherence", "snr"])} = ["all"]

        % opts.spatial_info {mustBeNumeric} = 0
        % opts.sparsity {mustBeNumeric} = 0
        % opts.mutual_info {mustBeNumeric} = 0
        % opts.entropy {mustBeNumeric} = 0
        % opts.kld {mustBeNumeric} = 0
        % opts.spatial_coherence {mustBeNumeric} = 0
        % opts.snr {mustBeNumeric} = 0
        % opts.all {mustBeNumeric} = 1
    end
    m = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% Prepare some data
    % convert to greyscale, needed for some metrics
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

%%%%%%%%%%%%%%%% Skaggs spatial information content (bits per second)
% Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
% https://onlinelibrary.wiley.com/doi/epdf/10.1002/%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K
% Markus et al. (1994) Spatial information content and reliability of hippocampal CA1 neurons: effects of visual input
% https://onlinelibrary.wiley.com/doi/epdf/10.1002/hipo.450040404
    if ismember("spatial_info",opts.metrics) || ismember("all",opts.metrics)
        pi = dmap ./ sum(dmap,'all','omitmissing'); % dwell time probability
        ro = sum(rmap(:) .* pi(:),'all','omitmissing'); % overall firing rate
        s = sum(pi(:) .* (rmap(:)./ro) .* log2(rmap(:)./ro),'all','omitmissing');
        m.skaggs_si_bits_per_sec = s;
    end

%%%%%%%%%%%%%%%% Skaggs spatial information content (bits per spike)
% From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
% The information rate given by formula (1) is measured in bits per second. If it is
% divided by the overall mean ring rate of the cell (expressed in spikes per second),
% then a different kind of information rate is obtained, in units of bits per spike|let us
% call it the information per spike. This is a measure of the specificity of the cell: the
% more grandmotherish" the cell, the more information per spike. 
    if ismember("spatial_info",opts.metrics) || ismember("all",opts.metrics)
        m.skaggs_si_bits_per_spike = s ./ ro;
    end

%%%%%%%%%%%%%%%% Sparsity
% From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
% The sparsity measure is  an adaptation  to space of  a formula invented by Treves and Rolls (1 99 1); the adaptation measures the 
% fraction of the environment  in which a cell  is active. Intuitively, a sparsity of, say, 0.1 means that the place field of the cell 
% occupies 1/10 of the area the rat traverses.
    if ismember("sparsity",opts.metrics) || ismember("all",opts.metrics)
        m.sparsity = ( sum(pi(:).*rmap(:),'all','omitmissing').^2 ) ./ ( sum(pi(:).*(rmap(:).^2),'all','omitmissing') );
    end

%%%%%%%%%%%%%%%% Mutual information
% Intuitively, mutual information measures the information that X and Y share: It measures how much knowing one of these variables reduces 
% uncertainty about the other. For example, if X and Y are independent, then knowing X does not give any information about Y and vice versa, 
% so their mutual information is zero. At the other extreme, if X is a deterministic function of Y and Y is a deterministic function of 
% X then all information conveyed by X is shared with Y: knowing X determines the value of Y and vice versa. As a result, in this case the 
% mutual information is the same as the uncertainty contained in Y (or X) alone, namely the entropy of Y (or X). Moreover, this mutual 
% information is the same as the entropy of X and as the entropy of Y. (A very special case of this is when X and Y are the same random variable.)
% https://en.wikipedia.org/wiki/Mutual_information
    if ismember("mutual_info",opts.metrics) || ismember("all",opts.metrics)
        m.mutual_info =  MutualInformation(rmap8(:),dmap8(:));
    end

%%%%%%%%%%%%%%%% Shannon's entropy (rmap only)
% help for matlab function entropy which does the same thing:
% E = entropy(I) returns E, a scalar value representing the entropy of an
% intensity image.  Entropy is a statistical measure of randomness that can be
% used to characterize the texture of the input image.  Entropy is defined as
% -sum(p.*log2(p)) where p contains the histogram counts returned from IMHIST.
    if ismember("entropy",opts.metrics) || ismember("all",opts.metrics)
        m.shannon_entropy = -sum(p.*log2(p));
    end

%%%%%%%%%%%%%%%% Cross entropy
% In information theory, the cross entropy between two probability distributions {\displaystyle p} p and {\displaystyle q} q over the same underlying 
% set of events measures the average number of bits needed to identify an event drawn from the set, if a coding scheme is used that is optimized for 
% an "unnatural" probability distribution {\displaystyle q} q, rather than the "true" distribution {\displaystyle p} p.
    if ismember("entropy",opts.metrics) || ismember("all",opts.metrics)
        ce0 = p .* log2(q); % using log base 2 means the output is measured in bits
        m.cross_entropy = -sum(ce0,'all','omitmissing');
    end

%%%%%%%%%%%%%%%% Kullback–Leibler divergence
% In mathematical statistics, the Kullback–Leibler divergence (also called relative entropy) is a measure of how one probability distribution diverges from 
% a second, expected probability distribution. In the simple case, a Kullback–Leibler divergence of 0 indicates that we can expect similar, if not the same, 
% behavior of two different distributions, while a Kullback–Leibler divergence of 1 indicates that the two distributions behave in such a different 
% manner that the expectation given the first distribution approaches zero. In simplified terms, it is a measure of surprise.
% https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
    if ismember("kld",opts.metrics) || ismember("all",opts.metrics)
        m.kldivergence = sum(p .* (log2(p) - log2(q)));
    end

%%%%%%%%%%%%%%%% Symmetric Kullback–Leibler divergence
% This quantity has sometimes been used for feature selection in classification problems, where P and Q are the conditional pdfs of a feature under two different classes.
    if ismember("kld",opts.metrics) || ismember("all",opts.metrics)
        KL1 = sum(p .* (log2(p) - log2(q)));
        KL2 = sum(q .* (log2(q) - log2(p)));
        m.kldivergence_symmetric = (KL1 + KL2) / 2;
    end

%%%%%%%%%%%%%%%% Jensen-Shannon divergence
% In probability theory and statistics, the Jensen–Shannon divergence is a method of measuring the similarity between two probability distributions.
% It is based on the Kullback–Leibler divergence, with some notable (and useful) differences, including that it is symmetric and it is always a finite value. 
% The square root of the Jensen–Shannon divergence is a metric often referred to as Jensen-Shannon distance.
% https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
    if ismember("kld",opts.metrics) || ismember("all",opts.metrics)
        logQvect = log2((q + p) / 2);
        jsd = 0.5 * (sum(p .* (log2(p) - logQvect)) + sum(q .* (log2(q) - logQvect)));
        m.jsdivergence = sqrt(jsd);
    end

%%%%%%%%%%%%%%%% Spatial Coherence 
% Cacucci et al. (2007)
% https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007
% The spatial coherence for each firing rate map was computed as the mean correlation between the firing rate of each bin with the 
% aggregate rate of the 24 nearest bins.
    if ismember("spatial_coherence",opts.metrics) || ismember("all",opts.metrics)
        meanf = ones([5 5]);
        meanf(3,3) = 0;
        U = imfilter(rmap,meanf,'replicate');
        m.spatial_coherence = corr(rmap(:),U(:),'rows','pairwise','type','Pearson');
    end

%% Signal to noise ratio / selectivity
% Called selectivity by Skaggs et al. (1996)
% The selectivity measure is  equal to the spatial maximum  fir- ing rate  divided  by  the mean  firing  rate  of the cell. The more tightly 
% concentrated the cell's activity, the higher  the selectivity. A cell  with  no spatial tuning at all will  have  a selectivity of  1; 
% there is  in principle  no upper  limit. A similar measure was  used by Barnes et al. (1983), except that the "out-of-field" firing rate was 
% used instead of the mean rate. The present definition is prefer- able because it does  not depend  on identifying  a  "place field," and because 
% it is much less sensitive to noise. 
    if ismember("snr",opts.metrics) || ismember("all",opts.metrics)
        m.signal_to_noise = max(rmap,[],"all",'omitmissing') ./ mean(rmap,'all','omitmissing');
    end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JointEntropy sub function
function H = JointEntropy(X)
    % JointEntropy: Returns joint entropy (in bits) of each column of 'X'
    % by Will Dwinnell
    %
    % H = JointEntropy(X)
    %
    % H = calculated joint entropy (in bits)
    % X = data to be analyzed
    %
    % Last modified: Aug-29-2006
    
    
    % Sort to get identical records together
    X = sortrows(X);
    
    % Find elemental differences from predecessors
    DeltaRow = (X(2:end,:) ~= X(1:end-1,:));
    
    % Summarize by record
    Delta = [1; any(DeltaRow')'];
    
    % Generate vector symbol indices
    VectorX = cumsum(Delta);
    
    % Calculate entropy the usual way on the vector symbols
    H = Entropy(VectorX);
    
    
    % God bless Claude Shannon.
    
    % EOF
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entropy sub function
function H = Entropy(X)
    % Entropy: Returns entropy (in bits) of each column of 'X'
    % by Will Dwinnell
    %
    % H = Entropy(X)
    %
    % H = row vector of calculated entropies (in bits)
    % X = data to be analyzed
    %
    % Example: Measure sample entropy of observations of variables with
    %   1, 2, 3 and 4 bits of entropy.
    %
    % Note: Estimated entropy values are slightly less than true, due to
    % finite sample size.
    %
    % X = ceil(repmat([2 4 8 16],[1e3,1]) .* rand(1e3,4));
    % Entropy(X)
    %
    % Last modified: Nov-12-2006
    
    % Establish size of data
    [n m] = size(X);
    
    % Housekeeping
    H = zeros(1,m);
    
    for Column = 1:m,
        % Assemble observed alphabet
        Alphabet = unique(X(:,Column));
	    
        % Housekeeping
        Frequency = zeros(size(Alphabet));
	    
        % Calculate sample frequencies
        for symbol = 1:length(Alphabet)
            Frequency(symbol) = sum(X(:,Column) == Alphabet(symbol));
        end
	    
        % Calculate sample class probabilities
        P = Frequency / sum(Frequency);
	    
        % Calculate entropy in bits
        % Note: floating point underflow is never an issue since we are
        %   dealing only with the observed alphabet
        H(Column) = -sum(P .* log2(P));
    end
    
    
    % God bless Claude Shannon.
    
    % EOF
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MutualInformation sub function
function I = MutualInformation(X,Y);
    % MutualInformation: returns mutual information (in bits) of the 'X' and 'Y'
    % by Will Dwinnell
    %
    % I = MutualInformation(X,Y);
    %
    % I  = calculated mutual information (in bits)
    % X  = variable(s) to be analyzed (column vector)
    % Y  = variable to be analyzed (column vector)
    %
    % Note: Multiple variables may be handled jointly as columns in matrix 'X'.
    % Note: Requires the 'Entropy' and 'JointEntropy' functions.
    %
    % Last modified: Nov-12-2006
    
    if (size(X,2) > 1)  % More than one predictor?
        % Axiom of information theory
        I = JointEntropy(X) + Entropy(Y) - JointEntropy([X Y]);
    else
        % Axiom of information theory
        I = Entropy(X) + Entropy(Y) - JointEntropy([X Y]);
    end
    
    
    % God bless Claude Shannon.
    
    % EOF
end
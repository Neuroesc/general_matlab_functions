function d = computeScatterDensity(v1,v2,varargin)
% computeScatterDensity for a set of XY points get their density
% For plotting scatter graphs colored by data density it can be helpful to have
% a vector which specifies the data density for each data point                                
%
% USAGE
%
% d = computeScatterDensity(v1,v2) process with default settings
%
% INPUT
%
% 'v1' - Mx1, dimension 1 coordinate data
%
% 'v2' - Mx1, dimension 2 coordinate data
%
% 'r1' - Scalar, resolution of density estimate in dimension 1
%        units are bins.
%
% 'r2' - Scalar, resolution of density estimate in dimension 2
%        units are bins.
%
% 'k' - Scalar, Gaussian smoothing SD of density estimate, set to 0 for no smoothing
%        units are bins.
%
% 'xlimits' - Mx2, [min max] bin edge limits in dimension 1
%        default is [lower bound - range/100, upper bound + range/100]
%
% 'ylimits' - Mx2, [min max] bin edge limits in dimension 2
%        default is [lower bound - range/100, upper bound + range/100]
%
% OUTPUT
%
% 'd' - Mx1, density of data in v1 and v2
%         units are proportion of data points
%
% NOTES
% 1. This is mainly for visualization, more specific and appropriate density
% measures will be needed for statistical analyses
%
% EXAMPLE
%
% figure
% tiledlayout(2,2)
% nexttile
% N = 1000; 
% v1 = normrnd(0,10,[N 1]);
% v2 = normrnd(0,40,[N 1]);
% d = computeScatterDensity(v1,v2);
% scatter(v1,v2,20,d,'filled','o','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
% axis square
% 
% nexttile
% N = 1000; 
% v1 = [normrnd(0,10,[round(N/2) 1]); normrnd(5,1,[N-round(N/2) 1])];
% v2 = normrnd(0,40,[N 1]);
% d = computeScatterDensity(v1,v2,'r1',64,'r2',64,'k',1);
% scatter(v1,v2,20,d,'filled','o','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
% axis square
% 
% nexttile
% N = 1000; 
% v1 = [normrnd(0,10,[round(N/2) 1]); normrnd(5,5,[N-round(N/2) 1])];
% v2 = [normrnd(0,10,[round(N/2) 1]); normrnd(5,5,[N-round(N/2) 1])];
% d = computeScatterDensity(v1,v2,'r1',128,'r2',128,'k',1);
% scatter(v1,v2,20,d,'filled','o','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
% axis square
% 
% nexttile
% N = 1000; 
% v1 = [normrnd(0,10,[round(N/3) 1]); normrnd(30,10,[N-round(N/3) 1])];
% v2 = [normrnd(0,10,[round(N/3) 1]); normrnd(30,10,[N-round(N/3) 1])];
% d = computeScatterDensity(v1,v2,'r1',64,'r2',64,'k',1);
% scatter(v1,v2,20,d,'filled','o','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
% axis square
% 
% SEE ALSO histcounts2, imgaussfilt, interp2

% HISTORY
%
% version 1.0.0, Release 04/04/25 Initial release, previous copy lost
% version 1.0.1, Release 05/04/25 Renamed, improved comments
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    p = inputParser;
    addRequired(p,'v1',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addRequired(p,'v2',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addParameter(p,'r1',64,@(x) isscalar(x) && isnumeric(x) && ~isnan(x));   
    addParameter(p,'r2',64,@(x) isscalar(x) && isnumeric(x) && ~isnan(x));   
    addParameter(p,'k',1,@(x) isscalar(x) && isnumeric(x) && ~isnan(x));  
    [a,b] = bounds(v1(:)); 
    c = range(v1(:));
    addParameter(p,'xlimits',[a-c/100,b+c/100],@(x) isnumeric(x) && ~any(isnan(x)));  
    [a,b] = bounds(v2(:)); 
    c = range(v2(:));
    addParameter(p,'ylimits',[a-c/100,b+c/100],@(x) isnumeric(x) && ~any(isnan(x)));  
    parse(p,v1,v2,varargin{:});
    config = p.Results;
    v1 = config.v1(:);
    v2 = config.v2(:);

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % bivariate histogram of data
    [f,ye,xe] = histcounts2(v2(:),v1(:),[config.r1 config.r2],'XBinLimits',config.ylimits,'YBinLimits',config.xlimits);
    if config.k>0
        f = imgaussfilt(f,config.k,'Padding',0,'FilterDomain','spatial');
    end
    f = f ./ sum(f,"all"); % normalize to unit sum
    
    % get bin centers
    xg = movmean(xe,2,'Endpoints','discard');
    yg = movmean(ye,2,'Endpoints','discard');
    [X,Y] = meshgrid(xg,yg);

    % interpolate to find density values for each data point
    d = interp2(X,Y,f,v1(:),v2(:),'linear');







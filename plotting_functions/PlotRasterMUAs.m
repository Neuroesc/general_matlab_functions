function [h1,h2,p] = PlotRasterMUAs(ax,spt,clu,varargin)
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
% n_clu = 12;
% n_spikes = 10;
% coef = 0.5;
% sigma = 1;
% spt = [];
% clu = [];
% for cc = 1:n_clu
%     spt = [spt; normrnd(cc*coef,sigma,n_spikes,1)];
%     clu = [clu; ones(n_spikes,1).*cc];
% end
%
% figure
% ax = gca;
% [h1,h2,p] = PlotRasterMUAs(ax,spt,clu,'marker_width',0.01,'marker_height',1,'marker_colors',cool(max(clu(:))),'linear_fit',1);
%
% SEE ALSO patch, polyfit

% HISTORY
%
% version 1.0.0, Release 05/01/26 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    % Parse inputs
    p = inputParser;   
    addRequired(p,'ax',@(x) ishandle(x) );      
    addRequired(p,'spt',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addRequired(p,'clu',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addParameter(p,'marker_width',0.01,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'marker_height',1,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'marker_colors',cool(max(clu(:))),@(x) isnumeric(x) );
    addParameter(p,'linear_fit',0,@(x) isscalar(x));       
    parse(p,ax,spt,clu,varargin{:});
    config = p.Results;

    if numel(config.clu) ~= numel(config.spt)
        error('N spikes (%d) does not equal N clusters (%d) exiting',numel(config.spt),numel(config.clu))
    end    
    if size(config.marker_colors,1) ~= numel(config.spt)
        warning('N spikes (%d) does not equal N colors (%d), using defaults',numel(config.spt),size(config.marker_colors,1))
        config.marker_colors = cool(max(config.clu(:)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % prepare vectors
    y = config.clu(:);
    x = config.spt(:);
    c = config.marker_colors(y,:);
    
    % remove NaNs
    nindx = isnan(x) | isnan(y);
    y = y(~nindx);
    x = x(~nindx);
    c = c(~nindx,:);
    
    % marker XY matrices
    X = [x'-config.marker_width; x'+config.marker_width; x'+config.marker_width; x'-config.marker_width];
    config.marker_height = config.marker_height ./ 2;
    Y = [y'-config.marker_height; y'-config.marker_height; y'+config.marker_height; y'+config.marker_height];  
        
    % plot data
    h1 = patch(ax,X,Y,'k'); 
    hold on;
    set(h,'facevertexcdata',c,'facecolor','flat','edgecolor','none') 
    
    % calculate linear fit line
    p = polyfit(x,y,1); % p(1) = slope, p(2) = intercept
    xr = linspace(ax.XLim(1),ax.XLim(2),100);
    yr = polyval(p,xr);
    
    axis manual
    if config.linear_fit
        h2 = plot(xr,yr,'k','LineWidth',1);
    end
    ax.YLim = [min(config.clu(:))-config.marker_height max(config.clu(:))+config.marker_height];















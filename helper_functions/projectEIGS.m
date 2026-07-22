function [Xs,Ys,Zs,F1] = projectEIGS(vecs,sres,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROJECTEIGS  projects unit vectors onto a unit sphere
% calculates the spherical KSDE of a population of unit vectors (place field eigenvectors or positional unit vectors)
% for every point on the unit sphere we calculate the spherical distance between every vector end point and the sphere
% point. These are then Gaussian weighted and summed. The final KDE is a map of these summed values.
%
% USAGE:
%         [Xs,Ys,Zs,F1] = projectEIGS(vecs,sres,sigma)
%
% INPUT:
%         vecs - unit vectors to be projected (Mx3)
%         sres - resolution (number of points) in spherical grid
%         sigma - the sigma/smoothing/SD of the Gaussian weighting, higher means a smoother output
%
% OUTPUT:
%    Xs,Ys,Zs - coordinates of the spherical grid, these are in a format that can be used directly by surf
%    F1 - the KDE map, values correspond to the coordinates in Sp
%
% EXAMPLES:
%
%     npoints = 1000; % number of random points
%     t = 2*pi*rand(npoints,1); % generate random azimuth
%     p = asin(-1+2*rand(npoints,1)); % generate random pitch
%     [X,Y,Z] = sph2cart(t,p,ones(npoints,1)); % convert to XYZ
%     [Xs,Ys,Zs,F1] = projectEIGS([X,Y,Z],64,10); % project these onto sphere
%     figure % plot the result
%     subplot(1,2,1)
%     plot3(X,Y,Z,'ko')
%     daspect([1 1 1])
%     axis xy tight
%     subplot(1,2,2)
%     surf(Xs,Ys,Zs,F1,'EdgeColor','none');
%     daspect([1 1 1])
%     axis xy tight
%
% See also: ACOSD DOT SPHERE

% HISTORY:
% version 1.0.0, Release 03/05/18 Initial release
% version 2.0.0, Release 14/11/18 Modified to use spherical gaussian
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
inps = {'sres','sigma'};
vals = {'64','10'};
reqd = [0 0];
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        if reqd(ff)
            error('ERROR: vargin %s missing... exiting',inps{ff});            
        end        
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
% get sphere coordinates
[Xs,Ys,Zs] = sphere(sres);

F1 = NaN(size(Xs));
if isempty(vecs)
    return
end

sig1 = (sqrt(2*pi) .* sigma);
sig2 = 1/sig1;

vpoints = vecs;
vpoints(any(isnan(vpoints),2),:) = [];
F1 = NaN(size(Xs));
for i = 1:numel(Xs)
    vsphere = [Xs(i) Ys(i) Zs(i)];
    vsphere = repmat(vsphere,length(vpoints(:,1)),1);
    dp = acosd(dot(vsphere,vpoints,2));
    dp_norm = (exp(-0.5 * (dp./sigma).^2) ./ sig1) ./ sig2;
    F1(i) = sum(dp_norm);
end
F1 = real(F1);

% figure
% subplot(1,2,1)
% plot3(eigs(:,1),eigs(:,2),eigs(:,3),'ko')
% daspect([1 1 1])
% axis xy tight
% 
% subplot(1,2,2)
% surf(Xs,Ys,Zs,F1,'EdgeColor','none');
% daspect([1 1 1])
% axis xy tight












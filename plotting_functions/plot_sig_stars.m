function [res,mres,sts,axb,h_b,h_t,h_bo,h_to] = plot_sig_stars(ds,gs,ax,opts)
% plot_sig_stars plot significance brackets on a graph
% Conducts statistical test on a data series, followed by multiple comparisons,
% then plots significance brackets (horizontal lines denoting significance) and
% asterisks to denote level of significance.
%
% USAGE
%
% plot_sig_stars(ax,ds,gs) process with default settings
%
% plot_sig_starsr(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'ax' - Axis handle, axis handle you would like to plot to
%       Default value is gca
%
% 'ds' - Numeric vector, data values provided as a vertical vector, these are
%       the values plotted in your graph
%
% 'gs' - Numeric vector, grouping values provided as a vertical vector, this is
%       typically the x-axis location data were plotted to, but can be any
%       value, data values will be grouped according to this vector and
%       significance brackets will be plotted according to their values. See
%       note 1.
%
% 'ds' - Numeric vector, data values provided as a vertical vector, these are
%       the values plotted in your graph
%
% 'test' - String, statistical test type that is desired, can be 'anova','kw'
%       for Kruskal-Wallis, or 'friedman' for Friedman's test, this is the
%       omnibus test.
%       Default value is 'anova'.
%
% 'display' - Logical scalar, set to true to see statistical results printed to
%       the command window
%       Default value is true.
%
% 'multcomp' - Logical scalar, set to true if you want significance brackets
%       plotted for post-hoc tests/multiple comparisons - this is what is
%       typically shown by significance brackets.
%       Default value is true.
%
% 'whiskers' - Logical scalar, set to true if you want the significance brackets
%       to have vertical 'whiskers' or caps at the ends, this is purely visual.
%       Default value is false.
%
% 'force_omnibus' - Logical scalar, set to true if you want the omnibus test
%       to be shown even if the test result is not significant (in this case the
%       p-value will be displayed as n.s.)
%       Default value is false.
%
% 'spacing_coef' - Numeric scalar, controls the vertical spacing between
%       significance brackets. A larger number means larger spacing and
%       corresponds to a percentage of the vertical distance spanned by the
%       y-axis.
%       Default value is 3
%
% 'horizon_coef' - Numeric scalar, controls the vertical spacing between
%       the bottom significance bracket and the top of the graph. A larger 
%       number means larger spacing and corresponds to a percentage of the 
%       vertical distance spanned by the y-axis.
%       Default value is 0, which means the first significance bracket will
%       overlap the top bar of the axis, if you want to use 'box on' this will
%       need to be increased.
%
% 'textshift' - Numeric scalar, controls the vertical spacing between
%       each significance bracket and the corresponding asterisks. A larger 
%       number means larger spacing and corresponds to a percentage of the 
%       vertical distance spanned by the y-axis.
%       Default value is 0, which means the text sits above the significance
%       bar, vertically centered, which for an asterisk plots it above the bar.
%
% 'whisker_coef' - Numeric scalar, controls the vertical span of any whiskers
%       plotted, if 'whiskers' is set to true. A larger number means larger 
%       whiskers and corresponds to a percentage of the vertical distance spanned 
%       by the y-axis.
%       Default value is 0.8
%
% 'font_size' - Numeric scalar, the font size of the asterisks and omnibus test
%       result plotted above significance brackets.
%       Default value is 10
%
% 'bracket_color' - Numeric vector, string, the color of the line used to draw the
%       significance bars, can be in Matlab color format, e.g. 'k', [0 0 0], 'r', 'red'
%       Default value is 'k', can be set to 'none' for invisible lines
%
% 'text_color' - Numeric vector, string, the color of the text used to give the
%       test results, can be in Matlab color format, e.g. 'k', [0 0 0], 'r', 'red'
%       Default value is 'k', can be set to 'none' for invisible text
%
% OUTPUT
%
% 'out' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% NOTES
% 1. The grouping vector must be the same length as the data vector, the numbers
%   in the grouping vector do not need to be integers, but they will need to be
%   consistent for values to be grouped together correctly. For example, a data
%   vector: 
%   [1 2 3 4 5 6 7 8 9 10] 
%   could have the grouping vector:
%   [1 1 1 1 1 2.2 2.2 2.2 2.2 2.2]
%   Two groups will be formed, one group for all the data values corresponding to
%   a group value of 1 and another for the values corresponding to a group value
%   of 2.2, furthermore, if there is a difference between these groups the
%   significance bars will be plotted between 1 and 2.2 on the x-axis, so these
%   values control both the grouping and the x-location of plotting.
%
% 2. An omnibus test is the test performed before multiple comparisons or
%   post-hoc tests, e.g. this function runs an ANOVA by default as the omnibus
%   test, only if this test is significant should multiple comparisons be
%   performed and this function will not run these tests if the omnibus test is
%   not significant.
%
% 3. In some instances you may want to plot the results of the omnibus test, for
%   example, if you have only 2 groups there will be no post-hoc tests because
%   they would just repeat the omnibus. In this case plot_sig_stars will plot a
%   significance bar using the output of the omnibus test instead.
% 
% This
%   behaviour can be disabled by setting 'omnibus' to false.
%
% EXAMPLE
%
% ds = reshape(rand(50,4).*[1:4],[],1);
% gs = reshape(ones(50,4).*[1:4],[],1);
% figure
% plot(gs,ds,'ko');
% box off
% ax = gca;
% [res,mres,sts,axb,h_b,h_t,h_bo,h_to] = plot_sig_stars(ds,gs,ax);
% 
% SEE ALSO anova1, multcompare

% HISTORY
%
% version 1.0.0, Release 15/04/26 Initial release
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
        ds {mustBeNumeric}
        gs {mustBeNumeric}
        ax matlab.graphics.axis.Axes = gca

        opts.test {mustBeMember(opts.test,{'anova','kw','friedman'})} = 'anova'
        opts.display logical = true
        opts.multcomp logical = true
        opts.whiskers logical = false
        opts.force_omnibus logical = false

        opts.spacing_coef {mustBeNumeric} = 2 % spacing between brackets, higher means more spacing
        opts.horizon_coef {mustBeNumeric} = 0 % spacing between bottom bracket and top of axis, higher means more spacing
        opts.textshift {mustBeNumeric} = 0 % spacing between text and brackets, higher means more spacing
        opts.whisker_coef {mustBeNumeric} = 0.8 % spacing between bottom bracket and top of axis, higher means more spacing

        opts.font_size {mustBeNumeric} = 10
        opts.bracket_color {mustBeA(opts.bracket_color, ["char", "string", "double"])} = 'k'
        opts.text_color {mustBeA(opts.text_color, ["char", "string", "double"])} = 'k'
    end
    axb = [];
    h_b = [];
    h_t = [];
    h_bo = [];
    h_to = [];
    mres = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% STATS
    % get some information
    ds = double( ds(:) );
    gs = double( gs(:) );
    [gu,~,gidx] = unique(gs);    
    ngroups = numel(gu);

    % group stats
    stats = struct;
    stats.meds = accumarray(gidx,ds,[],@(x) median(x,"all",'omitmissing'));
    stats.mens = accumarray(gidx,ds,[],@(x) mean(x,"all",'omitmissing'));
    stats.stds = accumarray(gidx,ds,[],@(x) std(x,[],"all",'omitmissing'));
    sem = @(x,d,m) std(x,d,'all',m) ./ sqrt( sum(~isnan(x),d,m) );
    stats.sems = accumarray(gidx,ds,[],@(x) sem(x,[],'omitmissing'));        
    nfun = @(x,d,m) sum(~isnan(x),d,m);
    stats.n = accumarray(gidx,ds,[],@(x) nfun(x,[],'omitmissing')); 

    display_text = sprintf('\t\tgroup n: %s\n\t\tmeans: %s\n\t\tmedians: %s\n\t\tstd: %s\n\t\tsem: %s',mat2str(stats.n'),mat2str(stats.mens(:)',3),mat2str(stats.meds(:)',3),mat2str(stats.stds(:)',3),mat2str(stats.sems(:)',3) );
    if opts.display
        disp( display_text )     
    end 

    %% Run omnibus test and comparisons
    % Omnibus is the overall test, such as the ANOVA that you run before running
    % multiple comparisons
    [res,sts,tomni] = local_stats_fun(ds,gs,opts);

    %% If more than 2 groups, run multiple comparisons
    % Thankfully multcompare can work out and adapt to the omnibus test that
    % was run
    p = [];
    if ngroups>2 && res{4}<=.05 % res{4} = omnibus p-value
        [p,gm,mres,tmulti] = local_multi_fun(sts,opts);
    end

%%%%%%%%%%%%%%%% SIG BRACKETS AND STARS
    % ordinarily we only want to plot post-hoc results if a series of criteria
    % are met, there are more than 2 groups (see below), the omnibus was
    % significant, at least one post-hoc test was significant, and the user
    % actually wants to plot post-hoc tests.
    print_multi = ngroups>2 && res{4}<=.05 && any(p<=.05) && opts.multcomp;

    % if there are only 2 groups and omnibus is significant
    % we want to plot the results of the omnibus because there will be no
    % post-hoc tests, but the groups differ, so we should plot a bracket    
    print_omni = false;
    if ngroups==2 && res{4}<=.05 
        print_omni = true;
    elseif opts.force_omnibus
        print_omni = true;
    end

    if print_multi || print_omni
        % create an invisible axis to plot the brackets above, this ensures
        % consistent spacing regardless of the axis values
        axb = axes('Position',ax.Position,'Color','none');
            axis(axb,'off')
            axb.YLim = [0 1];
            axb.XLim = ax.XLim;
            linkaxes([axb,ax],'x');
            hold(axb,'on')
            axis(axb,'manual')
    
        % work out bracket positions
        spacing = (opts.spacing_coef/100).*range(axb.YLim);
        horizon = axb.YLim(2) + ((opts.horizon_coef/100).*range(axb.YLim));
        tshift = (opts.textshift/100).*range(axb.YLim);

        if print_multi
            pidx = p<=.05;
            p = p(pidx);
            gm = gm(pidx,:);
            tmulti = tmulti(pidx);
        
            n_brackets = sum( p<=.05 );
            bracket_x_positions = gu(gm);
            bracket_y_positions = [horizon horizon] + spacing .* (0:(n_brackets-1))';

            if opts.whiskers
                bracket_x_positions = [bracket_x_positions(:,[1 1]), bracket_x_positions, bracket_x_positions(:,[end end])];
                bracket_y_positions = [bracket_y_positions(:,[1 1]), bracket_y_positions, bracket_y_positions(:,[end end])];
                v_drop = (opts.whisker_coef/100).*range(axb.YLim);
                bracket_y_positions = bracket_y_positions + [-v_drop +v_drop 0 0 -v_drop +v_drop];
                text_x = mean(bracket_x_positions(:,[3 4]),2,"omitmissing");
                text_y = mean(bracket_y_positions(:,[3 4]),2,"omitmissing");
            else
                text_x = mean(bracket_x_positions,2,"omitmissing");
                text_y = mean(bracket_y_positions,2,"omitmissing");
            end

            % plot brackets
            x_all = reshape( padarray(bracket_x_positions,[0,1],NaN,'post')', [], 1);
            y_all = reshape( padarray(bracket_y_positions,[0,1],NaN,'post')', [], 1);
            h_b = line(axb,x_all,y_all,'Color',opts.bracket_color,'Clipping','off');
        
            % add significance text
            h_t = text(axb,text_x,text_y+tshift,tmulti,'VerticalAlignment','middle','HorizontalAlignment','center','Clipping','off','FontSize',opts.font_size,'Interpreter','tex','Color',opts.text_color);

        elseif print_omni
            bracket_x_positions = [min(gu,[],'omitmissing') max(gu,[],'omitmissing')];
            bracket_y_positions = [horizon horizon] + spacing;

            if opts.whiskers
                bracket_x_positions = [bracket_x_positions(:,[1 1]), bracket_x_positions, bracket_x_positions(:,[end end])];
                bracket_y_positions = [bracket_y_positions(:,[1 1]), bracket_y_positions, bracket_y_positions(:,[end end])];
                v_drop = (opts.whisker_coef/100).*range(axb.YLim);
                bracket_y_positions = bracket_y_positions + [-v_drop +v_drop 0 0 -v_drop +v_drop];
                text_x = mean(bracket_x_positions(:,[3 4]),2,"omitmissing");
                text_y = mean(bracket_y_positions(:,[3 4]),2,"omitmissing");                
            else
                text_x = mean(bracket_x_positions,2,"omitmissing");
                text_y = mean(bracket_y_positions,2,"omitmissing");
            end

            % plot brackets
            x_all = reshape( padarray(bracket_x_positions,[0,1],NaN,'post')', [], 1);
            y_all = reshape( padarray(bracket_y_positions,[0,1],NaN,'post')', [], 1);
            h_bo = line(axb,x_all,y_all,'Color',opts.bracket_color,'Clipping','off');

            % add significance text
            % exponent = floor(log10(res{4}));
            % mantissa = res{4} / 10^exponent; 
            % t_omni = sprintf('{\\itp} = %.1f x 10^{%d}',mantissa,exponent);
            h_to = text(axb,text_x,text_y+tshift,tomni,'VerticalAlignment','bottom','HorizontalAlignment','center','Clipping','off','FontSize',opts.font_size,'Interpreter','tex','Color',opts.text_color);
        end
    end
    
    axes(ax); % return original axis as current

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OMNIBUS STATS FUNCTION
function [res,sts,t] = local_stats_fun(ds,gs,opts)
    switch opts.test
        case {'anova'}
            [~,a,sts] = anova1(ds,gs,'off');           
            eta = a{2,2} ./ a{4,2}; % Eta-squared = SSbetween / SStotal
            res = {a{2,3} a{3,3} a{2,5} a{2,6} eta};
            exponent = floor(log10(res{4}));
            mantissa = res{4} / 10^exponent;         
            if res{4}<=.05
                stest = '(significant)';
            else
                stest = '(not significant)';
            end            
            display_text = sprintf('\t\tANOVA: F(%d,%d) = %.1f, p = %.1f x 10^%d, eta2 = %.4f %s',res{1},res{2},res{3},mantissa,exponent,res{5},stest);
                           
        case {'kw'}
            [~,a,sts] = kruskalwallis(ds,gs,'off');
            epsilon = a{2,5} ./ ( (sum(~isnan(ds(:))).^2-1)/(sum(~isnan(ds(:)))+1) ); % epsilon squared
            res = {a{2,3} a{3,3} a{2,5} a{2,6} epsilon};  
            exponent = floor(log10(res{4}));
            mantissa = res{4} / 10^exponent; 
            if res{4}<=.05
                stest = '(significant)';
            else
                stest = '(not significant)';
            end              
            display_text = sprintf('\t\tKruskal-wallis: X2(%d,%d) = %.1f, p = %.1f x 10^%d, epsilon2 = %.4f %s',res{1},res{2},res{3},mantissa,exponent,res{5},stest );
 
        case {'friedman'}
            [~,gidx] = sort(gs,'ascend');
            ds2 = ds(gidx);
            ds3 = reshape(ds2,[],ngroups);
            [~,a,sts] = friedman(ds3,1,'off');              
            res = {a{2,3} a{3,3} a{2,5} a{2,6}};   
            exponent = floor(log10(res{4}));
            mantissa = res{4} / 10^exponent;
            if res{4}<=.05
                stest = '(significant)';
            else
                stest = '(not significant)';
            end              
            display_text = sprintf('\t\tFriedman: X2(%d,%d) = %.1f, p = %.1f x 10^%d %s',res{1},res{2},res{3},mantissa,exponent,stest);

        % case {'repeatt'}
        %     [~,gidx] = sort(gs,'ascend');
        %     ds2 = ds(gidx);
        %     ds3 = reshape(ds2,[],ngroups);
        %     [~,omnip,~,sts] = ttest(ds3(:,1),ds3(:,2));
        %     meds = accumarray(gidxd,ds,[],@nanmean);      
        %     stds = accumarray(gidxd,ds,[],@nanstd);  
        %     sems = accumarray(gidxd,ds,[],@nansem);  
        %     res = [sts.df sts.tstat omnip meds(:)'];      
        %     if opts.display
        %         disp(sprintf('\t\tPaired ttest: t(%d) = %.1f, p = %s\n\t\tgroup n: %s\n\t\tmeans: %s\n\t\tSD: %s\n\t\tSEM: %s',sts.df,sts.tstat,mat2str(omnip,3),mat2str([sum(~isnan(ds3(:,1))) sum(~isnan(ds3(:,2)))]),mat2str(meds(:)',3),mat2str(stds(:)',3),mat2str(sems(:)',3) ) )
        %         if opts.exactp
        %             disp(omnip);
        %         end            
        %     end    
    end      

    t = 'n.s.';
    if res{4}<=.05
        t = '*';
        if res{4}<=.01
            t = '**';                    
            if res{4}<=.001
                t = '***';   
            end
        end
    end

    % display results
    if opts.display
        disp( display_text )     
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MULTCOMPARE STATS FUNCTION
function [p,g,c,t] = local_multi_fun(sts,opts)
    [c,~,~,~] = multcompare(sts,'display','off','ctype','dunn-sidak');
    p = c(:,end);
    g = c(:,1:2);
    x = sort(unique(g));

    if opts.display
        for pp = 1:length(p)
            if p<=.05
                stest = '(significant)';
            else
                stest = '(not significant)';
            end
            disp( sprintf('\t\t\tgroup %d v %d: p = %.2e %s',x(g(pp,1)),x(g(pp,2)),p(pp),stest) ) 
        end
    end 

    t = cell(length(p),1);
    for i = 1:length(p)
        if p(i)<=.05
            t(i) = {'*'};
            if p(i)<=.01
                t(i) = {'**'};                    
                if p(i)<=.001
                    t(i) = {'***'};   
                end
            end
        end
    end
end















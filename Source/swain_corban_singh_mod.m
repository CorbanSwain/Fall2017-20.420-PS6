function swain_corban_singh_mod

end

%% Corban Swain Utilities

function [val, S] = qgetfield(S, fieldname, default)
if ~isfield(S, fieldname)
    S.(fieldname) = default;
end
val = S.(fieldname);
end

function fighand = makefigure(fs)
    function val = qgf(fieldname, default)
        [val, fs] = qgetfield(fs, fieldname, default);
    end
fighand = setupfig(fs.n, fs.title, qgf('position', 'auto'));
qgf('sub', [1, 1]);
if prod(fs.sub) < length(fs.plots)
    error('The subplot dimensions are too small for the number of plots.');
end
for i = 1:length(fs.plots)
    subplot(fs.sub(1), fs.sub(2), i);
    hold on;
    makeplot(fs.plots{i});
end
end

classdef Plotter
    properties
        X
        Y
        LineSpec
        ColorOrderLen
        CycleLen
        MarkerSize
        MarkerFaceColor
        LineWidth
        Grid = 'off'
        XLabel = 'x'
        YLabel = 'y'
        XScale = 'auto'
        YScale = 'auto'
        XLim = 'auto'
        YLim = 'auto'
        Legend
        LegendTitle
        LegendLineWidth
        LegendLoc
    end
    methods
        function makeplot(ps)
            function val = qgf(fieldname, default)
                [val, ps] = qgetfield(ps, fieldname, default);
            end
            isf = @(fieldname) isfield(ps, fieldname);
            if ~all(isfield(ps, {'x', 'y'}))
                error('X and Y values must be passed in the plotspec struct.');
            end
            if iscell(ps.x) && iscell(ps.y)
                n = length(ps.x);
                cycle = 1;
                qgf('linespec', {'-'});
                if  n ~= length(ps.y)
                    error('x and y cell arrays have mismatched dimensions.');
                end
                for i = 1:n
                    fieldname = 'color_order_len';
                    if isf(fieldname)
                        if (mod(i, ps.(fieldname)) == 1) && (i ~= 1)
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                        end
                    end
                    fieldname = 'cycle_len';
                    if isf(fieldname)
                        if (mod((i - 1), ps.(fieldname)) == 0) && (i ~= 1)
                            cycle = cycle + 1;
                        end
                    end
                    % FIXME - Implement linespec like markersize
                    % FIXME - maybe functionalize this more? Be able to take in shorter
                    % cell arrays then vars.
                    p = plot(ps.x{i}, ps.y{i}, ps.linespec{cycle});
                    fieldname = 'markersize';
                    if isf(fieldname)
                        sizespec = ps.(fieldname){i};
                        if ~isempty(sizespec)
                            p.MarkerSize = sizespec;
                        end
                    end
                    fieldname = 'markerfacecolor';
                    if isf(fieldname)
                        colorspec = ps.(fieldname){i};
                        if ~isempty(colorspec)
                            p.MarkerFaceColor = colorspec;
                        end
                    end
                    
                    fieldname = 'linewidth';
                    if isf(fieldname)
                        lw = ps.(fieldname){i};
                        if ~isempty(lw)
                            p.LineWidth = lw;
                        end
                    end
                    
                end
            else
                if iscell(ps.x) || iscell(ps.y)
                    error('x and y values must both be matrices or cell arrays.');
                end
                % FIXME - Check for Linespec as cell array ... actually use single
                % value above under all conditions if a cell array is not passed. ...
                % maybe convert x and y to a 1 X 1 cell array to keep everything in the
                % same place...
                plot(ps.x, ps.y, qgf('linespec', '-'));
            end
            ax = p.Parent;
            grid(ax, qgf('grid', 'off'));
            ax.XScale = qgf('xscale', 'linear');
            ax.YScale = qgf('yscale', 'linear');
            xlabel(qgf('xlabel','x'));
            ylabel(qgf('ylabel','y'));
            xlim(qgf('xlim','auto'));
            ylim(qgf('ylim','auto'));
            fieldname = 'legend';
            if isf(fieldname')
                leg = legend(ps.(fieldname));
                fieldname = 'legend_title';
                if isf(fieldname')
                    title(leg, ps.(fieldname));
                end
                leg.LineWidth = qgf('legend_linewidth', 0.5);
                leg.Location = qgf('legend_loc','best');
            end
        end
    end
end

function new_fig = setupfig(n,title,location)
% SETUPFIGURE Sets up a new figure.
%
% n: figure number
% title: figure title
% location: figure position, [left bottom width height]
%
new_fig = figure(n); clf; hold on;
box off; grid on;
new_fig.Name = title;
if nargin > 2
    if location ~= 'auto'
        new_fig.Position = location;
    end
end
end

function savefig(fig,fig_name)
% SAVEFIGURE Saves the passed figure as a 300 dpi png.

if ~isdir([pwd filesep 'Figures'])
    mkdir 'Figures'
end
f = gobjects(1,1);
name = '';
switch nargin
    case 0
        f = gcf;
    case 1
        f = fig;
    case 2
        f = fig;
        name = fig_name;
end
if isempty(name)
    if isempty(f.Name)
        name = 'Untitled';
    else
        name = fig.Name;
    end
else
    if ~isempty(f.Name)
        name = [name, '-', f.Name];
    end
end
filename = ['Figures' filesep name];
print(f,filename,'-dpng','-r300');
end

function corban_figure_defaults
% CORBANFIGUREDEFAULTS Sets default values to make pretty figures.
fontSize = 13;
font = 'Helvetica';
set(groot, ...
    'defaultLineMarkerSize', 40,...
    'defaultLineLineWidth', 3, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultAxesTitleFontWeight', 'normal', ...
    'defaultAxesFontName', font, ...
    'defaultAxesLabelFontSizeMultiplier', 1.1, ...
    'defaultAxesLineWidth', 2, ...
    'defaultFigureColor', [1 1 1], ...
    'defaultTextInterpreter', 'tex', ...
    'defaultTextFontSize',fontSize, ...
    'defaultTextFontName', font ...
    );
end

function cleanup
clc;
clear;
end
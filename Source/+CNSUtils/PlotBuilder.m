classdef PlotBuilder
    properties
        X
        Y
        LineSpecCycleLength
        ColorOrderCycleLength
        LineSpec = {'-'} % setter update cell length
        MarkerSize % setter update cell length
        MarkerFaceColor % setter update cell length
        LineWidth % setter update cell length
        Grid = 'off'
        XScale = 'linear'
        YScale = 'Linear'
        XLabel = 'x'
        YLabel = 'y'
        XLim = 'auto'
        YLim = 'auto'
        LegendLabels
        LegendLineWidth = 1.5
        LegendLocation = 'best'
        LegendTitle
    end
    methods
        function obj = set.X(obj, val)
            obj.X = convert2cell(val);
        end
        function obj = set.Y(obj, val)
            obj.Y = convert2cell(val);
        end
        function obj = PlotBuilder(X, Y)
            if nargin == 2
                obj.X = X;
                obj.Y = Y;
            end
        end
        
        function plot(obj)
            if isempty(obj.Y)
                error('Y-Values must be passed to build a plot.');
            end
            nYVals = length(obj.Y);
            if isempty(obj.X)
                for i = 1:nYVals
                    obj.X{i} = 1:length(obj.Y{i});
                end
            end
            if  nYVals ~= length(obj.X)
                error('x and y cell arrays have mismatched dimensions.');
            end
            plotSettingNames = {'MarkerSize','MarkerFaceColor', ...
                                'LineWidth'};
            lineSpecIndex = 1;
            for i = 1:nYVals
                % FIXME - maybe functionalize this more? Be able to take in
                % shorter cell arrays then vars.
                p = plot(obj.X{i}, obj.Y{i}, obj.LineSpec{lineSpecIndex});
                if i == 1, axes = p.Parent; end
                
                if ~isempty(obj.ColorOrderCycleLength)
                    if (mod(i, obj.ColorOrderCycleLength) == 0)
                        axes.ColorOrderIndex = 1;
                    end
                end
                
                if ~isempty(obj.LineSpecCycleLength)
                    if (mod(i, obj.LineSpecCycleLength) == 0)
                        lineSpecIndex = lineSpecIndex + 1;
                    end
                end
                
                for iSetting = 1:length(plotSettingNames)
                    name = plotSettingNames{iSetting};
                    propertyVal = obj.(name);
                    if ~isempty(propertyVal) && ~isempty(propertyVal{i})
                        p.(name) = propertyVal;
                    end
                end
           
            end
        grid(axes, obj.Grid);
        axes.XScale = obj.XScale;
        axes.YScale = obj.YScale;
        xlabel(obj.XLabel);
        ylabel(obj.YLabel);
        xlim(obj.XLim);
        ylim(obj.YLim);
        if ~isempty(obj.LegendLabels)
            legendHandle = legend(obj.LegendLabels);
            if ~isempty(obj.LegendTitle)
                title(legendHandle, obj.LegendTitle);
            end
            legendHandle.LineWidth = obj.LineWidth;
            legendHandle.Location = obj.LegendLocation;
        end
    end
end
end
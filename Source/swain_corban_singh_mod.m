function swain_corban_singh_mod



%% Global Variableas
cacheName = "swain_SimulationCacheMod";


% vv CHANGE THIS to 10 to recreate figures shown in report
simTime = 0.5; % s


atpConc = 2e-3; % M
dt0 = 2e-4; % s
axonCoverage = 0.6; % percentage
mtSpacing = 20; % nm
dyneinPerMt = 2;
nVesicles = 100;
        

%% Functions

    function simulateAndSave
        [dT, dX, dY, vT, vX, vY, vB] = mod_simulate(simTime, ...
            atpConc, dt0, axonCoverage, mtSpacing, dyneinPerMt, nVesicles);
        [dT2, dX2, dY2, vT2, vX2, vY2, vB2] = mod_simulate(simTime,...
            atpConc, dt0, axonCoverage, mtSpacing, 0, nVesicles);
        save(cacheName, 'dT', 'dX', 'dY', 'vT', 'vX', 'vY', 'vB', ...
            'dT2', 'dX2', 'dY2', 'vT2', 'vX2', 'vY2', 'vB2');
    end


%% Figures
figures{11} = @testfig1;
    function fb = testfig1
        nFig = 11;
        dt = 1e-4;
        nSteps = 1e5;
        pos = zeros(nSteps, 2); 
        for iStep = 2:nSteps
            prevPos = pos(iStep - 1, :); 
            pos(iStep, :) = rand_walk(prevPos, ModConstants.D, dt, ...
                [0 Inf], [0 2e3]);
        end
        pb = CNSUtils.PlotBuilder;
        pb.X{1} = pos(:, 1);
        pb.Y{1} = pos(:, 2);
        pb.PointColor{1} = (1:nSteps) .* dt;
        pb.LineWidth{1} = 0.5; 
        pb.Marker{1} = 'none';
        pb.MarkerFaceColor{1} = 'flat';        
        fb = CNSUtils.FigureBuilder;
        fb.Number = nFig;
        fb.PlotBuilders = pb;
    end

figures{12} = @testfig2;
    function fb = testfig2
        nFig = 12;
        dt = 1E-2;
        nSteps = 1e5;
        nTrials = 30;
        finalPos = zeros(nTrials, 2);
        initialPos = zeros(1, 2);
        parfor iTrial = 1:nTrials            
            prevPos = initialPos;
            for iStep = 2:nSteps
                prevPos = rand_walk(prevPos, ModConstants.D, dt);
            end
            finalPos(iTrial, :) = prevPos;
        end
        pb = CNSUtils.PlotBuilder;
        for iTrial = 1:nTrials            
            pb.X{iTrial} = [initialPos(1), finalPos(iTrial, 1)];
            pb.Y{iTrial} = [initialPos(2), finalPos(iTrial, 2)];
        end
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = nFig;
        fb.PlotBuilders = pb;
    end

figures{13} = @testFigure3;
    function fb = testFigure3
        nFig = 7;
        load(cacheName, 'vX', 'vX2', 'vY', 'vY2', 'vT', 'vT2');
        pb = CNSUtils.PlotBuilder;
        maxX1 = max(max(vX));
        maxX2 = max(max(vX2));
        maxX = max([maxX1, maxX2]);
        for iVes = 1:nVesicles
            pb.X{iVes} = vX(:, iVes) ./ 1e3;
            pb.Y{iVes} = (vY(:, iVes) + 2000 .* (iVes)) ./ 1e3;
            pb.PointColor{iVes} = vT(:, iVes); 
            pb.ColorMap{iVes} = 'jet';
%             pb.CAxisLimits{iVes} = [9 10];
            pb.LineWidth{iVes} = 1;
        end
        for iVes = 1:nVesicles
            plotInd = iVes + nVesicles;
            pb.X{plotInd} = vX2(:, iVes) ./ 1e3;
            pb.Y{plotInd} = (vY2(:, iVes) + 2000 .* (plotInd + 10)) ...
                ./ 1e3;
            pb.PointColor{plotInd} = vT2(:, iVes); 
            pb.ColorMap{plotInd} = 'jet';
%             pb.CAxisLimits{plotInd} = [9 10];
            pb.LineWidth{plotInd} = 1;
        end
        pb.ColorBar = 'on';
        pb.XLabel = 'x-position (\mum)';
        pb.YLabel = 'Vessicle Number';
        pb.YTicks = [1 202 222 423];
        pb.YTickLabels = {'1', '100', '1', '100'}; 
        pb.EdgeColorMethod = 'flat';
        pb.XLim = [0 maxX] ./ 1e3;
        pb.YLim = [0 (2000 * (((nVesicles + 11)* 2) ))] ./ 1e3;
        fb = CNSUtils.FigureBuilder;
        fb.Number = nFig;
        fb.PlotBuilders = pb;
        fb.Position = [245 83 1403 736];
    end

figures{7} = @modFigure1;
    function fb = modFigure1
        fignum = 7;
        pb = CNSUtils.PlotBuilder;
        pb.X = 1:110;
        pb.Y = ModConstants.K_ON_EFF(1:110);
        pb.YScale = 'log';
        pb.XLabel = 'Vessicle-Dynein Distance (nn)';
        pb.YLabel = 'K_{on} (M^{-1} s^{-1})';
        pb.YLim = [min(pb.Y{1}), max(pb.Y{1})];
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf('%d - Binding Affinity vs Distance', fignum);
        fb.Position = [70 493 888 414];
        fb.PlotBuilders = pb;
    end

figures{8} = @modFigure2;
    function fb = modFigure2
        halfWidth = ModConstants.AXON_WIDTH ./ 2 .* axonCoverage;
        dInitialPos(:, 2) = colon(-halfWidth, mtSpacing, halfWidth) ...
            + (ModConstants.AXON_WIDTH ./ 2);
        nMt = size(dInitialPos, 1);
        nFig = 8;
        load(cacheName, 'vX', 'vX2', 'vY', 'vY2', 'vT', 'vT2');
        pb = CNSUtils.PlotBuilder;
        maxX1 = max(max(vX));
        maxX2 = max(max(vX2));
        maxX = max([maxX1, maxX2]);
        for i = 1:nMt
            pb.X{i} = [0 maxX] ./ 1e3;
            pb.Y{i} = [dInitialPos(i, 2), dInitialPos(i, 2)] ./ 1e3 + 0.15;
            pb.LineSpec{i} = 'k-';
            pb.LineWidth{i} = 0.5;
        end
        fudge = 0.0;
        axonLims = [0-fudge 2-fudge 2.3-fudge 4.3+fudge] + 0.15;
        for i = 1:length(axonLims)
            n = nMt + i;
            pb.X{n} = [0 maxX] ./ 1e3;
            pb.Y{n} = [axonLims(i) axonLims(i)];
            pb.LineSpec{n} = 'k:';
            pb.LineWidth{n} = 2;
        end
        for iVes = nMt + length(axonLims) + 1
            id = 40;
            pb.X{iVes} = vX(:, id) ./ 1e3;
            pb.Y{iVes} = (vY(:, id)) ./ 1e3 + 0.15;
            pb.PointColor{iVes} = vT(:, id); 
            pb.ColorMap{iVes} = 'jet';
%             pb.CAxisLimits{iVes} = [9 10];
            pb.LineWidth{iVes} = 2.5;
        end
        for iVes = nMt + length(axonLims) + 2
            id = 31;
            plotInd = iVes;
            pb.X{plotInd} = vX2(:, id) ./ 1e3;
            pb.Y{plotInd} = (vY2(:, id) + 2300) ./ 1e3 + 0.15;
            pb.PointColor{plotInd} = vT2(:, id); 
            pb.ColorMap{plotInd} = 'jet';
%             pb.CAxisLimits{plotInd} = [9 10];
            pb.LineWidth{plotInd} = 2.5;
        end
        
        pb.LineSpecCycleLength = 1;
        pb.ColorBar = 'on';
        pb.XLabel = 'x-position (\mum)';
        pb.YLabel = 'y-position (\mum)';
        pb.YTicks = [0 2 2.3 4.3] + 0.15;
        pb.YTickLabels = {'0', '2', '0', '2'}; 
        pb.EdgeColorMethod = 'flat';
        pb.XLim = [0 maxX] ./ 1e3;
        pb.YLim = [0 4550] ./ 1e3;
        fb = CNSUtils.FigureBuilder;
        fb.Name = sprintf(strcat("%d - Single Vessicle Movement ", ...
            "with and wo Dynein"), nFig);
        fb.Number = nFig;
        fb.PlotBuilders = pb;
        fb.Position = [245 83 1403 736];
    end

figures{9} = @modFigure3;
    function fb = modFigure3
        nFig = 9;
        load(cacheName, 'vX', 'vX2', 'vY', 'vY2', 'vT', 'vT2');
        pb = CNSUtils.PlotBuilder;
        maxX1 = max(max(vX));
        maxX2 = max(max(vX2));
        maxX = max([maxX1, maxX2]);
        
        fudge = 0.01;
        axonLims = [0-fudge 2-fudge 2.3-fudge 4.3+fudge] + 0.15;
        for i = 1:length(axonLims)
            n = i;
            pb.X{n} = [0 maxX] ./ 1e3;
            pb.Y{n} = [axonLims(i) axonLims(i)];
            pb.LineSpec{n} = 'k:';
            pb.LineWidth{n} = 2;
        end
        
        for iVes = 1:nVesicles
            n = iVes + length(axonLims);
            pb.X{n} = vX(:, iVes) ./ 1e3;
            pb.Y{n} = (vY(:, iVes)) ./ 1e3 + 0.15;
            pb.PointColor{n} = vT(:, iVes); 
            pb.ColorMap{n} = 'jet';
%             pb.CAxisLimits{iVes} = [9 10];
            pb.LineWidth{n} = 1;
        end
        for iVes = 1:nVesicles
            plotInd = iVes + nVesicles + length(axonLims);
            pb.X{plotInd} = vX2(:, iVes) ./ 1e3;
            pb.Y{plotInd} = (vY2(:, iVes) + 2300) ./ 1e3 + 0.15;
            pb.PointColor{plotInd} = vT2(:, iVes); 
            pb.ColorMap{plotInd} = 'jet';
%             pb.CAxisLimits{plotInd} = [9 10];
            pb.LineWidth{plotInd} = 1;
        end
        pb.ColorBar = 'on';
        pb.XLabel = 'x-position (\mum)';
        pb.YLabel = 'y-position (\mum)';
        pb.YTicks = axonLims;
        pb.YTickLabels = {'0', '2', '0', '2'}; 
        pb.EdgeColorMethod = 'flat';
        pb.XLim = [0 maxX] ./ 1e3;
        fb = CNSUtils.FigureBuilder;
        fb.Name = sprintf("%d - Vessicle Movement with and wo Dynein", ...
            nFig);
        fb.Number = nFig;
        fb.PlotBuilders = pb;
        fb.Position = [348 34 1403 736];
    end


%% Main Block

    function main
        CNSUtils.cleanup;
        
        simulateAndSave;
        
        CNSUtils.FigureBuilder.setDefaults;
        for i = 7:9
            figure(figures{i}()); drawnow;
            CNSUtils.FigureBuilder.saveFigure;
        end
    end % function main

tic; main; toc;
end

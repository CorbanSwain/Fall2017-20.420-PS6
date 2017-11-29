function swain_corban_singh_mod



%% Global Variableas



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

figures{7} = @modFigure1;
    function fb = modFigure1
        nFig = 7;
        simTime = 10; % s
        atpConc = 2e-3; % M
        dt0 = 2e-4; % s
        axonCoverage = 0.6; % percentage
        mtSpacing = 20; % nm
        dyneinPerMt = 2;
        nVesicles = 100;
        cacheName = "SimulationCache-HighD.mat";
        doLoad = true;
        if doLoad
            load(cacheName);
        else
            [dT, dX, dY, vT, vX, vY, vB] = mod_simulate(simTime, atpConc, ...
                dt0, ...
                axonCoverage, mtSpacing, dyneinPerMt, nVesicles);
            [dT2, dX2, dY2, vT2, vX2, vY2, vB2] = mod_simulate(simTime,...
                atpConc, dt0, ...
                axonCoverage, mtSpacing, 0, nVesicles);
            save(cacheName, 'dT', 'dX', 'dY', 'vT', 'vX', 'vY', 'vB', ...
                'dT2', 'dX2', 'dY2', 'vT2', 'vX2', 'vY2', 'vB2');
        end
        
        pb = CNSUtils.PlotBuilder;
        maxX1 = max(max(vX));
        maxX2 = max(max(vX2));
        maxX = max([maxX1, maxX2]);
        for iVes = 1:nVesicles
            pb.X{iVes} = vX(:, iVes);
            pb.Y{iVes} = vY(:, iVes) + 2000 .* (iVes);
            pb.PointColor{iVes} = vB(:, iVes); 
            pb.ColorMap{iVes} = 'cool';
            pb.CAxisLimits{iVes} = [0 1];
            pb.LineWidth{iVes} = 1;
        end
        for iVes = 1:nVesicles
            plotInd = iVes + nVesicles;
            pb.X{plotInd} = vX2(:, iVes);
            pb.Y{plotInd} = vY2(:, iVes) + 2000 .* (plotInd + 10);
            pb.PointColor{plotInd} = vB2(:, iVes); 
            pb.ColorMap{plotInd} = 'cool';
            pb.CAxisLimits{plotInd} = [0 1];
            pb.LineWidth{plotInd} = 1;
        end
        pb.EdgeColorMethod = 'flat';
        pb.XLim = [0 maxX];
        pb.YLim = [0 (2000 * ((nVesicles * 2) + 11))];
        fb = CNSUtils.FigureBuilder;
        fb.Name = "Vessicle Movement with Dynein (bound)";
        fb.Number = nFig;
        fb.PlotBuilders = pb;
    end

%% Main Block

    function main
%         CNSUtils.cleanup;
        CNSUtils.FigureBuilder.setDefaults;
        figure(figures{7}()); drawnow;
        CNSUtils.FigureBuilder.saveFigure;
        commandwindow;
    end % function main

tic; main; toc;
end

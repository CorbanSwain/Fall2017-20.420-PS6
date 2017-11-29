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
        simTime = 5; % s
        atpConc = 1e-3; % M
        dt0 = 1e-4; % s
        axonCoverage = 0.5; % percentage
        mtSpacing = 100; % nm
        dyneinPerMt = 1;
        nVesicles = 10;
        [dT, dX, dY, vT, vX, vY] = mod_simulate(simTime, atpConc, dt0, ... 
            axonCoverage, mtSpacing, dyneinPerMt, nVesicles);
        fb = CNSUtils.FigureBuilder;
    end

%% Main Block

    function main
        CNSUtils.cleanup;
        CNSUtils.FigureBuilder.setDefaults;
        figures{7}();
        drawnow; commandwindow;
    end % function main

tic; main; toc;
end

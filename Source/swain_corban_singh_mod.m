function swain_corban_singh_mod



%% Global Variableas



%% Figures
figures{11} = @testfig1;
    function fb = testfig1
        dt = 1e-3;
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
        fb.Number = 1;
        fb.PlotBuilders = pb;
    end

figures{12} = @testfig2;
    function fb = testfig2
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
        fb.Number = 2;
        fb.PlotBuilders = pb;
    end


%% Main Block

    function main
        CNSUtils.FigureBuilder.setDefaults;
        figure(figures{11}());
        drawnow; commandwindow;
    end % function main

tic; main; toc;
end

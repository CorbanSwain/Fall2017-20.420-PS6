function swain_corban_singh_mod



%% Global Variableas



%% Figures
figures{11} = @testfig1;
    function fb = testfig1
        dt = 1E-2;
        nSteps = 1e6;
%         pos = zeros(2, nSteps);
        
        for iStep = 2:nSteps
%             prevPos = pos(:, iStep - 1); 
%             pos(:,iStep) = rand_walk(prevPos, ModConstants.D, dt);
            prevPos = rand_walk(prevPos, ModConstants.D, dt);
        end
        pb = CNSUtils.PlotBuilder;
        pb.X = pos(1,[1, end]);
        pb.Y = pos(2,[1, end]);
        
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
        figure(figures{12}());
    end % function main

tic; main; toc;
end

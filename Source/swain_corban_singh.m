function swain_corban_singh
% SWAIN_CORBAN_SINGH Implementation of Singh, 2005 dynein model.



%% Global Variables
STALL_TIME = 1; % s
SIM_TIME = 4; % s



%% Figures

figures{2} = @fig2;
    function fb = fig2
        fignum = 2;
        atpConcs = [1e-3, 400e-6];
        nTrials = length(atpConcs);
        [T, X] = simulate(nTrials, atpConcs, ...
                          @SinghConstants.restoringForce, ...
                          STALL_TIME, SIM_TIME);
        
        pb = CNSUtils.PlotBuilder;
        pb.X = T;
        pb.Y = X;
        
        pb.X{2} = pb.X{2} + 0.5;
        pb.YLabel = 'Position (nm)';
        pb.XLabel = 'Time(s)';
        pb.XLim = [0 4];
        pb.LegendLabels = {'1 mM ATP', '400 \muM ATP'};
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf('%d - Typical Dynein Simulations', fignum);
        fb.PlotBuilders = pb;   
    end % function fig2

figures{3} = @fig3;
    function fb = fig3
        fignum = 3;
        atpConcs = [(1:10) .* 100e-6, (1.5:0.5:4) .* 1e-3]; % M
        force = @SinghConstants.restoringForce;
        nTrials = length(atpConcs);
        nRepeats = 50;
        stallForce = zeros(nRepeats, nTrials);
        fprintf('\tBeginning Simulation Loop, %d Repeats ...\n', ...
                nRepeats);
        parfor iRepeat = 1:nRepeats
%             Dynein.calcCache
            fprintf('%2d.', iRepeat);
            [~, X] = simulate(nTrials, atpConcs, force, ...
                             STALL_TIME, SIM_TIME);
            for iTrial = 1:nTrials
                stallForce(iRepeat,iTrial) = force(X{iTrial}(end));
            end % loop through trials
        end % repeat experiment loop
        fprintf('\n\tDone!\n');
        pb = CNSUtils.PlotBuilder;
        pb.X{1} = atpConcs .* 1e6; % uM
        pb.Y{1} = mean(stallForce);
        pb.YError{1} = std(stallForce, 0, 1);
        pb.XLabel = 'ATP Concentration (\muM)';
        pb.YLabel = 'Stall Force (pN)';
        pb.LineSpec = {'o-'};
        pb.MarkerSize = {6};
        pb.MarkerFaceColor = 'w';
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf('%d - Stalling Force vs ATP Concentration', ...
                          fignum);
        fb.PlotBuilders = pb;
    end % function fig3



%% Main Block

    function main
        % CNSUtils.cleanup;
        % close all;
        CNSUtils.FigureBuilder.setDefaults;
        figsToRun = 3;
        Dynein.calcCache;
        for iFig = figsToRun
            fprintf('\nRunning Figure %d\n', iFig);
            fb = figures{iFig}();
            fb = figure(fb);
            save(fb); 
        end
    end
tic
main;
toc
end


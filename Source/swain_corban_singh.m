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
        fprintf('\tBeginnning Simulations ... ')
        [T, X] = simulate(nTrials, atpConcs, ...
                          @SinghConstants.restoringForce, ...
                          STALL_TIME, SIM_TIME);
        fprintf('Done!\n');
        
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
        fb.Position = [5 664 568 300];
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
        fprintf('\tBeginning Simulation Loop, %d Repeats ... \n', ...
                nRepeats);
        parfor iRepeat = 1:nRepeats
%             Dynein.calcCache
            fprintf('%3d.', iRepeat);
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
        pb.MarkerSize = {8};
        pb.MarkerFaceColor = {'w'};
        pb.LineWidth = {2.5};
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf('%d - Stalling Force vs ATP Concentration', ...
                          fignum);
        fb.Position = [5 241 834 334];
        fb.PlotBuilders = pb;
    end % function fig3

figures{5} = @fig5;
    function fig5
        fignum = 5;
        loads = 0:0.05:1;
        atpConcs = [100e-6, 2e-3];
        nLoads = length(loads);
        nAtpConcs = length(atpConcs);
        nRepeats = 10;
        velocities = zeros(nAtpConcs, nLoads, nRepeats);
        simTime = 10; % s
        stallTime = 5;
        fprintf('\tBeginning Simulation Loop, %d Repeats ...\n', ...
                nRepeats);
        velocSlice = zeros(nLoads,nRepeats);
        for iAtpConc = 1:nAtpConcs
            fprintf('\tATP Concentration %d of %d ... \n', ...
                iAtpConc, nAtpConcs);
            atp = atpConcs(iAtpConc);
            parfor iRepeat = 1:nRepeats
                fprintf('%3d.', iRepeat);    
                [T, X] = simulate(nLoads, atp, loads, ...
                                  stallTime, simTime);
                vss = zeros(nLoads, 1);
                for iLoad = 1:nLoads
                    vss(iLoad) = X{iLoad}(end) ...
                        ./ T{iLoad}(end);
                end
                velocSlice(:,iRepeat) = vss;
                
            end
            fprintf('\n');
            velocities(iAtpConc,:,:) = shiftdim(velocSlice, -1);
        end % main for loop
        fprintf('\tDone!\n');
        
        pb = CNSUtils.PlotBuilder;
        n = 1;
        pb.X{n} = loads;
        pb.Y{n} = mean(velocities(n, :, :), 3);
        pb.YError{n} = std(velocities(n, :, :), 0, 3);
        pb.YLabel{n} = 'Velocity (nm/s), 100\muM ATP';
        pb.YLim{n} = [0, 180];
        n = n + 1;
        pb.X{n} = loads;
        pb.Y{n} = mean(velocities(n, :, :), 3);
        pb.YError{n} = std(velocities(n, :, :), 0, 3);
        pb.YLabel{n} = 'Velocity (nm/s), 2 mM ATP';
        pb.YLim{n} = [0 650];
        pb.AxisAssignment = [1, 2];
        pb.XLabel = 'Load (pN)';                
        pb.LineSpec = {'o-'};
        pb.MarkerFaceColor = {'w','w'};
        pb.MarkerSize = {8, 8};
        pb.LegendLabels = {'100 \muM ATP', '2 mM ATP'};
        pb.LineWidth = {2.5, 2.5};
        pb.Box = 'on';
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf('%d - Average Velocity vs. Load', fignum);
        fb.Position = [843 569 640 392];
        fb.PlotBuilders = pb;
    end

figures{6} = @fig6;
    function fb = fig6
        fignum = 6;
        loads = {0, 0.2, 0.4, 0.6}; % pN
        c = [1, 2, 5] .* 1e-6; % M
        atpConcs{1} = [c, ...
                       c * 10, ...
                       c * 100, ...
                       c * 1e3];
        atpConcs{2} = atpConcs{1};
        atpConcs{3} = atpConcs{1}(5:end);
        atpConcs{4} = atpConcs{1}(6:end);
        
        nLoads = length(loads);
        nRepeats = 10;
        
        velocities = cell(1, nLoads);        
        simTime = 20; % s
        stallTime = 10; %s
        fprintf('\tBeginning Simulation Loop, %d Repeats ...\n', ...
                nRepeats);
        for iLoad = 1:nLoads
            fprintf('\tATP Concentration %d of %d ... \n', ...
                iLoad, nLoads);
            load = loads{iLoad};
            atp = atpConcs{iLoad};
            nSims = length(atp);
            vSlice = zeros(nRepeats, nSims);
            parfor iRepeat = 1:nRepeats
                fprintf('%3d.', iRepeat);    
                [T, X] = simulate(nSims, atp, load, ...
                                  stallTime, simTime);
                v = zeros(1, nSims);
                for iSim = 1:nSims
                    v(iSim) = X{iSim}(end) ./ T{iSim}(end);
                end
                vSlice(iRepeat, :) = v;
            end
            fprintf('\n');
            velocities{iLoad} = vSlice;
        end % main for loop
        fprintf('\tDone!\n');
        
        pb = CNSUtils.PlotBuilder;
        for iLoad = 1:nLoads
            pb.X{iLoad} = atpConcs{iLoad} .* 1e6;
            pb.Y{iLoad} = mean(velocities{iLoad}, 1);
            pb.YError{iLoad} = std(velocities{iLoad}, 0, 1);
        end
        pb.YLabel = 'Velocity (nm/s)';
        pb.YLim = [1e-2, 1e3];
        pb.YScale = 'log';
        pb.XLabel = 'ATP Concentration (pN)';  
        pb.XLim = [0.5, 1e4];
        pb.XScale = 'log';
        pb.LineSpecCycleLength = 1;
        pb.LineSpec = {'o-', 's-', '^-', 'p-'};
        pb.MarkerFaceColor = {'w', 'w', 'w', 'w'};
        pb.MarkerSize = {8, 8, 8, 8};
        pb.LegendLabels = {'None', '0.2 pN', '0.4 pN', '0.6 pN'};
        pb.LegendTitle = 'Load';
        pb.LineWidth = {2.5, 2.5, 2.5, 2.5};
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf(['%d - Average Velocity vs. ', ...
                           'ATP Concentration'], ...
                          fignum);
        fb.Position = [];
        fb.PlotBuilders = pb;
    end % function fig 6



%% Main Block

    function main
        CNSUtils.cleanup;
        fprintf('Beginning Script.\n');
        
        if isempty(gcp('nocreate'))
            parpool('local');
        end
        
        CNSUtils.FigureBuilder.setDefaults;
        figsToRun = [6];
        for iFig = figsToRun
            fprintf('\nRunning Figure %d\n', iFig);
            fb = figures{iFig}();
            fb = figure(fb);
            save(fb); 
        end % figure for loop
        fprintf('\nScript Complete!\n\n');
    end % main

tic; main; toc;
end


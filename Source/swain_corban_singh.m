function swain_corban_singh
% SWAIN_CORBAN_SINGH Implementation of Singh, 2005 dynein model.



%% Global Variables
stallTime = 1; % s
simTime = 4; % s



%% Figures
figures{2} = @fig2;
    function fb = fig2
        fignum = 2;
        atpConcs = [1e-3, 400e-6];
        nTrials = length(atpConcs);
        
        T = zeros(1e5, nTrials);
        X = T;
        dyneins(1,nTrials) = Dynein;
        
        continueSim = true(1, nTrials);
        i = 1;
        % FIXME - this is very unclear
        while any(continueSim)
            simXs = X(:,continueSim);
            simTs = T(:,continueSim);
            simDyneins = dyneins(continueSim);
            
            simForces = SinghConstants.restoringForce(simXs(i,:));
            simAtpConcs = atpConcs(continueSim);
            
            [simDyneins, dt] = simDyneins.update(simAtpConcs, simForces);
            
            i = i + 1;
            simTs(i,:) = dt + simTs((i - 1), :);
            simXs(i,:) = [simDyneins.Position];
            simContinueSim = shouldContinueSim(simTs,simXs, i);
            
%             fprintf('Force %13.5f\n', simForces(1));
%             disp(simDyneins(1));
            T(:,continueSim) = simTs;
            X(:,continueSim) = simXs;
            dyneins(continueSim) = simDyneins;
            continueSim(continueSim) = simContinueSim;
        end % while simulation loop
        
        pb = CNSUtils.PlotBuilder;
        pb.X = cell(1, nTrials);
        pb.Y = pb.X;
        for iTrial = 1:nTrials
            lengthSim = find(T(2:end, iTrial) == 0, 1);
            pb.X{iTrial} = T(1:lengthSim, iTrial);
            pb.Y{iTrial} = X(1:lengthSim, iTrial);
        end
        pb.YLabel = 'Position (nm)';
        pb.XLabel = 'Time(s)';
        pb.LegendLabels = {'1 mM ATP', '400 \muM ATP'};
        
        fb = CNSUtils.FigureBuilder;
        fb.Number = fignum;
        fb.Name = sprintf('Typical Dynein Simulations');
        fb.PlotBuilders = pb;   
    end % function fig2



%% Functions
    function val = shouldContinueSim(T, X, last)
        T = T(1:last,:);
        X = X(1:last,:);
        nTrials = size(T, 2);
        val = true(1, nTrials);
        for i = 1:nTrials
            t = T(:,i);
            if t(end) > simTime
                val(i) = false;
            else
                x = X(:,i);
                stallTimePoint = t(end) - stallTime;
                if stallTimePoint > 0
                    iCheck = find(t > stallTimePoint, 1);
                    xCheck = x(iCheck);
                    if ~isempty(xCheck)
                        if abs(xCheck - x(end)) < eps
                            val(i) = false;
                        end
                    end % checking empty
                end % stall time point
            end % checking simTime
        end % for
    end % function stopCriterion(...)



%% Main
    function main
        CNSUtils.cleanup;
        % close all;
        CNSUtils.FigureBuilder.setDefaults;
        figsToRun = 2;
        for iFig = figsToRun
            fb = figures{iFig}();
            fb = figure(fb);
            save(fb);
        end
    end

main;
end


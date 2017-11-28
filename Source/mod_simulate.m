function [T, X] = mod_simulate(nTrials, atpConcs, force, simTime)



%% Input Handling
if isscalar(atpConcs)
    atpConcs = ones(1, nTrials) .* atpConcs;
else
    if numel(atpConcs) ~= nTrials
        error(['atpConcs must be a scalar or '
            'matrix with numel = nTrials.']);
    end
end




%% Setup
preAllocLength = 7e2 .* simTime;
[tempT, tempX] = deal(zeros(preAllocLength, nTrials));
dyneins(1,nTrials) = Dynein;



%% Main Simulation Loop
iStep = 1;
% FIXME - this is very unclear
while any(continueSim)
    % CHECK FOR FORCES
    simForces = force(continueSim);
    
    [simDyneins, dt] = simDyneins.update(atpConcs, forceFunc);
    
    iStep = iStep + 1;
    simTs(iStep,:) = dt + simTs((iStep - 1), :);
    simXs(iStep,:) = [simDyneins.Position];
    simContinueSim = shouldContinueSim(simTs,simXs, iStep, ...
        stallTime, simTime);

    tempT(:,continueSim) = simTs;
    tempX(:,continueSim) = simXs;
    dyneins(continueSim) = simDyneins;
    continueSim(continueSim) = simContinueSim;
end % while simulation loop



%% Trimming Simulation Timecourses
T = cell(1, nTrials);
X = T;
for iTrial = 1:nTrials
    lengthSim = find(tempT(2:end, iTrial) == 0, 1);
    T{iTrial} = tempT(1:lengthSim, iTrial);
    X{iTrial} = tempX(1:lengthSim, iTrial);
end
end % function simulate



function val = shouldContinueSim(T, X, last, stallTime, simTime)
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
end % function shouldContinueSim(...)
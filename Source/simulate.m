function [T, X] = simulate(nTrials, atpConcs, force, ...
    stallTime, simTime)



%% Input Handling
if isscalar(atpConcs)
    atpConcs = ones(1, nTrials) .* atpConcs;
else
    if numel(atpConcs) ~= nTrials
        error(['atpConcs must be a scalar or '
            'matrix with numel = nTrials.']);
    end
end

if isa(force, 'function_handle')
    forceFunction = force;
    forceOpt = 1;
elseif isscalar(force)
    forceOpt = 2;
else
    if numel(force) ~= nTrials
        error(['force must be a function handle, scalar, or '
            'matrix with numel = nTrials.']);
    end
    forceOpt = 3;
    force = force(:);
end



%% Setup
preAllocLength = 1e3 .* simTime;
tempT = zeros(preAllocLength, nTrials);
tempX = tempT;
dyneins(1,nTrials) = Dynein;
continueSim = true(1, nTrials);



%% Main Simulation Loop
iStep = 1;
% FIXME - this is very unclear
while any(continueSim)
    simXs = tempX(:,continueSim);
    simTs = tempT(:,continueSim);
    simDyneins = dyneins(continueSim);
    
    switch forceOpt
        case 1
            simForces = forceFunction(simXs(iStep,:));
        case 2
            simForces = force;
        case 3
            simForces = force(continueSim);
    end
    
    simAtpConcs = atpConcs(continueSim);
    
    [simDyneins, dt] = simDyneins.update(simAtpConcs, simForces);
    
    iStep = iStep + 1;
    simTs(iStep,:) = dt + simTs((iStep - 1), :);
    simXs(iStep,:) = [simDyneins.Position];
    simContinueSim = shouldContinueSim(simTs,simXs, iStep, ...
        stallTime, simTime);
    
    %             fprintf('Force %13.5f\n', simForces(1));
    %             disp(simDyneins(1));
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
end

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
end % function stopCriterion(...)
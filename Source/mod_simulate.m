function [dT, dX, dY, vT, vX, vY] = mod_simulate(simTime, atpConc, dt0, ...
    axonCoverage, mtSpacing, dyneinPerMt, nVesicle)
% axonCoverage: a percentage between 0 and 1


%% Input Handling


%% Dyenin Setup
preAllocLength = 7e2 .* simTime;
halfWidth = ModConstants.AXON_WIDTH ./ 2 .* axonCoverage;
dInitialPos(:, 2) = colon(-halfWidth, mtSpacing, halfWidth);
nMt = size(dInitialPos, 1);
nDynein = nMt * dyneinPerMt;
dInitialPos = [dInitialPos; zeros(nDynein - nMt, 2)];
if dyneinPerMt > 1
    for iDyn = 2:dyneinPerMt
        rows = (1:nMt) + (nMt .* (iDyn - 1));
        x0 = ModConstants.L_DYN .* (iDyn - 1);
        dInitialPos(rows, 1) = x0;
        dInitialPos(rows, 2) = dInitialPos(1:nMt, 2);
    end
end
dyneins = ModDynein(dInitialPos);
dyneinUpdateCount = ones(1, nDynein);
[tempDT, tempDX] = deal(zeros(preAllocLength, nDynein)); 



%% Vesicle Setup
timeVector = 0:dt0:simTime;
nTimePoints = length(timeVector);
vT(:, nVesicle) = timeVector;
[vX, vY] = deal(zeros(size(vT)));
vInitialPos(:,2) = linspace(0, ModConstants.AXON_WIDTH, nVesicle);
vesicles = Vesicle(vInitialPos);



%% Main Simulation Loop
for iTime = 2:nTimePoints
    time = timeVector(iTime);
    fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    % Dynein Update and Cache
    for iDyn = 1:nDynein
        
        while dyneins(iDyn).Time < time
%             fprintf("Updating Dynein %d\n", iDyn);
            dyneinUpdateCount(iDyn) = dyneinUpdateCount(iDyn) + 1;
            upCount = dyneinUpdateCount(iDyn);
            dyneins(iDyn) = dyneins(iDyn).update(atpConc);
            tempDT(upCount, iDyn) = dyneins(iDyn).Time;
            tempDX(upCount, iDyn) ... 
                = dyneins(iDyn).Position(1);
        end
    end
    
    % Vessicle Update and Cache
    vesicles = vesicles.update(dt0, time, dyneins);
    for iVes = 1:nVesicle
        pos = vesicles(iVes).Position;
        vX(iTime, iVes) = pos(1);
        vY(iTime, iVes) = pos(2);
    end
    fprintf(" Time: %8.5f", time);
end



%% Trimming Simulation Timecourses
[dT, dX, dY] = deal(cell(1, nDynein));
for iDyn = 1:nDynein
    uCount = dyneinUpdateCount(iDyn);
    dT{iDyn} = tempT(1:uCount, iDyn);
    dX{iDyn} = tempX(1:uCount, iDyn);
    dY{iDyn} = ones(uCount, 1) .* vInitialPos(iDyn, 2);
end
end % function simulate

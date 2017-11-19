function dt = mcTimeStep(rates)
%MCTIMESTEP Monte carlo time step.
%
%MCTIMESTEP(rates) Takes a vector or rates and randomly picks a time step b
%ased on those rates. Useful for Monte Carlo simulations.
if ~isvector(rates)
    error('rates must be a vector.')
end
dt = random('exp', 1 / sum(rates));
end


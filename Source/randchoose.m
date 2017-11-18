function choice = randchoose(weights, nTimes)
switch nargin
    case 1
        nTimes = 1;
    case 2
    otherwise
        warning('Too many arguments passed to randchoose.');
end
probs = weights ./ sum(weights, 2);
cumSumProbs = cumsum(probs, 2);
cumSumProbs = repmat(cumSumProbs, nTimes, 1);
nWeights = size(weights, 1);
nChoices = size(cumSumProbs, 1);
randomVar = rand(1, nChoices);
choice = zeros(nWeights, nTimes);
for iChoice = 1:nChoices
    r =  randomVar(iChoice);
    cs = cumSumProbs(iChoice, :);
    choice(iChoice) = find(r < cs, 1);
end
end
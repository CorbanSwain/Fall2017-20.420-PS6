function choice = randchoose(weights, nTimes)
switch nargin
    case 1
        nTimes = 1;
    case 2
    otherwise
        error('Unexpected number of arguments.');
end
binary = false;
if length(weights) == 1
    if weights > 1 || weights < 0
        error('A single weight must be a probability, 0 <= p <= 1.')
    end
    binary = true;
    probs = [weights, (1 - weights)];
else
    if any(isinf(weights))
        infTerms = (weights == Inf);
        weights(:) = 0;
        weights(infTerms) = 1;
    end
    wSum = sum(weights);
    probs = weights / wSum;
end
nWeightSets = size(weights, 1);
nChoices = nWeightSets * nTimes;
randomVar = rand(1, nChoices);
choice = zeros(nWeightSets, nTimes);
iProb = 0;
for iChoice = 1:nChoices
    iProb = iProb + 1;
    iProb(iProb > nWeightSets) = 1;
    r =  randomVar(iChoice);
    counter = 0;
    singleChoice = 0;
    cs = 0;
    while ~singleChoice
        counter = counter + 1;
        cs = cs + probs(iProb, counter);
        if r < cs, singleChoice = counter; end 
    end
    choice(iChoice) = singleChoice;
end
if binary
    choice = (choice == 1);
end
end
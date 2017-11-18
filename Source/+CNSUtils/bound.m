function val = bound(val, min, max, name)
%CHECKNUMBER Summary of this function goes here
%   Detailed explanation goes here
switch nargin
    case 3
        warn = false;
    case 4
        warn = true;
    otherwise
        error('Unexpected number of arguments.');
end
lessThanMin = val < min;
greaterThanMax = val > max;

if warn
    if all(~lessThanMin)
        warning('%s cannot be smaller than %d', name, min);
    end
    if all(~greaterThanMax)
        warning('%s cannot be larger than %d', name, max);
    end
end

val(lessThanMin) = min;
val(greaterThanMax) = max;
end


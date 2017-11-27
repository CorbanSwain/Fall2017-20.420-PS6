function pos = rand_walk(pos, D, dt, xlim, ylim)
angle = rand(1) .* 2 .* pi;
d = sqrt(D .* 4 .* dt);
dPos = [cos(angle), sin(angle)] .* d;
pos = pos + dPos;
switch nargin
    case 3

    case 5
        lim = [xlim', ylim'];
        while true
            limDiff = lim - pos;
            if all((limDiff .* [-1; 1]) >= 0), break; end
            if limDiff(1, 1) > 0
                dPos(1) = limDiff(1, 1) * 2;
            elseif limDiff(2, 1) < 0
                dPos(1) = limDiff(2, 1) * 2;
            end
            if limDiff(1, 2) > 0
                dPos(2) = limDiff(1, 2) * 2;
            elseif limDiff(2, 2) < 0
                dPos(2) = limDiff(2, 2) * 2;
            end
            pos = pos + dPos;
        end 
    otherwise
        % FIXME -  This should be at the beginning of the function...
        error('Unexpected number of arguments');
end
end
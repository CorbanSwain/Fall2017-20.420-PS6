function pos = rand_walk(pos, D, dt)
angle = rand(1) .* 2 * pi;
d = sqrt(D .* 4 .* dt);
pos(1) = pos + cos(angle) .* d;
pos(2) = pos + sin(angle) .* d;
end
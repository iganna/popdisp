function segments = noise_segments(segments, scale)

[c, ia, ib] = unique([segments(:,1:2); segments(:,3:4)], 'rows')
np = size(segments, 1)
points = [ib(1:np), ib((np+1):(2*np))]


% scale = 0.03
for i = 1:length(ia)
    delta = (randn(1,2)-0.5) * scale;
    s = [segments(points(:,1) == i, 1:2); segments(points(:,2) == i, 3:4)];
    s = s(1,:) + delta;
    segments(points(:,1) == i, 1:2) = repmat(s, sum(points(:,1) == i), 1);
    segments(points(:,2) == i, 3:4) = repmat(s, sum(points(:,2) == i), 1);
end
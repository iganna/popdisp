function points = segment_points(segments)
[c, ia, ib] = unique([segments(:,1:2); segments(:,3:4)], 'rows')
np = size(segments, 1)
points = [ib(1:np), ib((np+1):(2*np))]
function segments = split_segments(segments, dx)

[c, ia, ib] = unique([segments(:,1:2); segments(:,3:4)], 'rows')
np = size(segments, 1)
points = [ib(1:np), ib((np+1):(2*np))]

% dx = 0.1
n_seg = size(segments, 1)
idx_seg_del = [];
for i = 1:n_seg
    seg_d = norm(segments(i,3:4) - segments(i,1:2));
    if floor(seg_d / dx) < 2
        continue
    end
    n_split = floor(seg_d / dx);
    seg_d = segments(i,3:4) - segments(i,1:2);
    for j = 1:n_split
        if j == n_split
            seg_tmp = [segments(i,1:2) + seg_d/n_split * (j-1), segments(i,3:4)];
        else
            seg_tmp = [segments(i,1:2) + seg_d/n_split * (j-1), segments(i,1:2) + seg_d/n_split*j];
        end
        
        segments = [segments;seg_tmp];
    end
    idx_seg_del = [idx_seg_del; i];
end
segments(idx_seg_del,:) = [];


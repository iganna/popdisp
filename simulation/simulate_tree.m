
path_sim = 'sim_data/'
n_sim = 100

for i_sim = 1:n_sim

    path_tmp = [path_sim 's' num2str(i_sim), '_']
    % mkdir(path_tmp)
    %%

    n_leaves = 30
    x = rand(n_leaves)
    d = x*x'
    d = d - diag(diag(d))

    t = seqneighjoin(d)

    t_loc = phytree(get(t, 'Pointers'),  abs(randn(2*n_leaves - 1,1)))
    plot(t_loc, 'Type', 'equaldaylight')

%     saveas(gcf,[path_tmp 'tree_init.pdf'])


    %% Segments
    h = findall(gca, 'Type', 'Line')
    idx_branches = [];
    segments = [];
    for i = 1:length(h)
        if length(h(i).XData) == 2
            idx_branches = [idx_branches; i];
            segments = [segments;[h(i).XData(1) h(i).YData(1) h(i).XData(2) h(i).YData(2)]];
        end
    end
    y = h(1).YData;
    x = h(1).XData;

    close all  % close figure

    % Plot Segments
    f1 = figure; hold on;
    for j = 1:size(segments, 1)
        plot([segments(j,1) segments(j,3)], [segments(j,2)  segments(j,4)], '-')
    end
    axis equal
    saveas(gcf,[path_tmp 'segments_1_init.pdf'])
    close(f1)

    % Split segments
    segments = split_segments(segments, 0.1)

    f1 = figure; hold on;
    for j = 1:size(segments, 1)
        plot([segments(j,1) segments(j,3)], [segments(j,2)  segments(j,4)], '-')
    end
    axis equal
    saveas(gcf,[path_tmp 'segments_2_split.pdf'])
    close(f1)

    % Remember segments
    segments_init = segments;
%     writematrix(segments, [path_tmp 'segments_2_split.txt'])


    %% Point + noize

    segments = segments_init

    scale = 0.03
    segments = noise_segments(segments, scale)

    f1 = figure; hold on;
    for j = 1:size(segments, 1)
        plot([segments(j,1) segments(j,3)], [segments(j,2)  segments(j,4)], '-')
    end
    axis equal
    saveas(gcf,[path_tmp 'segments_3_noise.pdf'])
    close(f1)


%     writematrix(segments, [path_tmp 'segments_3_noise.txt'])

    %% moving by gradient

    xmin = min([segments(:,1); segments(:,3)]);
    xmax = max([segments(:,1); segments(:,3)]);
    ymin = min([segments(:,2); segments(:,4)]);
    ymax = max([segments(:,2); segments(:,4)]);
    step = 0.1
    [X, Y] = meshgrid((xmin):step:(xmax), (ymin):step:(ymax));


    n_peaks = 4
    m = [X(1,randi(size(X, 2), n_peaks,1))', Y(randi(size(Y, 1), n_peaks,1),1)]
    m = [m; [mean([segments(:,1); segments(:,3)]),mean([segments(:,2); segments(:,4)])]];
    n_peaks = size(m, 1);

    writematrix(m, [path_tmp 'peaks.txt'])

    syms x y
    n_peaks = 5
    z = 0;

    gscale = 1
    for i = 1:n_peaks
        z = z + exp(-(x-m(i,1))^2 / gscale-(y-m(i,2)).^2 / gscale);
    end
    g = gradient(z);


%     step = 0.1
%     [X, Y] = meshgrid(xmin:step:xmax, ymin:step:ymax);
    Z = double(subs(z, {x, y}, {X, Y}));   
%     f1 = figure; 
%     contour(X,Y,Z, 50);
%     axis equal
%     saveas(gcf,[path_tmp 'geo_1.pdf'])
%     close(f1)


    nstep = 10
    gstep = 0.2;
    segments_new = segments;
    points = segment_points(segments_new);
    for istep = 1:nstep    
        istep
        for i = 1:size(points, 1)
            pos = [segments_new(points(:,1) == i, 1:2); segments_new(points(:,2) == i, 3:4)];
            pos = pos(1,:);
            delta = (-double(subs(g, {x, y}, {pos(1), pos(2)})) * gstep)';
            pos = pos + delta;
            segments_new(points(:,1) == i, 1:2) = repmat(pos, sum(points(:,1) == i), 1);
            segments_new(points(:,2) == i, 3:4) = repmat(pos, sum(points(:,2) == i), 1);
        end
        segments_new = split_segments(segments_new, 0.1);
        points2 = segment_points(segments_new);

        segments_new = noise_segments(segments_new, 0.01);
        points = segment_points(segments_new);

        points = segment_points(segments_new);
    end
    writematrix(segments_new, [path_tmp 'segments_4_opt.txt'])
    writematrix(points, [path_tmp 'points_4_opt.txt'])

%     % Plot new tree on geography
%     f1 = figure; hold on;
%     contour(X,Y,Z, 50);
%     % for j = 1:size(segments, 1)
%     %     plot([segments(j,1) segments(j,3)], [segments(j,2)  segments(j,4)], 'k-')
%     % end
% 
%     for j = 1:size(segments_new, 1)
%         plot([segments_new(j,1) segments_new(j,3)], [segments_new(j,2)  segments_new(j,4)], 'r-')
%     end
%     axis equal
%     saveas(gcf, [path_tmp 'geo_2_tree_opt.pdf'])
%     close(f1)

    
    f1 = figure; hold on;
    for j = 1:size(segments_new, 1)
        plot([segments_new(j,1) segments_new(j,3)], [segments_new(j,2)  segments_new(j,4)], 'k-')
    end
    axis equal
    saveas(gcf, [path_tmp 'segments_4_opt.pdf'])
    close(f1)



    %% get points and their coordinates

    id_leaves = [];
    len_route = zeros(max(points(:)), 1);
    p_from = points(1,1);
    p_last = [];
    while length(p_from) > 0
        p_from_current = p_from(1);
        idx_to = [points(:,1) == p_from_current];
        idx_to(1) = false;

        d = ((segments_new(idx_to,1) - segments_new(idx_to,3)).^2 + ...
            (segments_new(idx_to,2) - segments_new(idx_to,4)).^2) .^ (1/2);
        if(sum(d>1)>0)
            'stop'
            break
        end
        sum(len_route(p_from_current,:))
        len_route(points(idx_to,2),:) = repmat(len_route(p_from_current,:), sum(idx_to), 1);

        len_route = [len_route, zeros(max(points(:)), 1)];
        len_route(points(idx_to,2), size(len_route, 2)) = d;



        p_last = [p_last; p_from(1)];
        p_from(1) = [];
        p_from = [p_from; points(idx_to,2)];

        if length(intersect(p_from, p_last)) > 0
            'aaa'
            break
        end
        if length(points(idx_to,2)) == 0
            id_leaves = [id_leaves, p_from_current];
        end
    end

    sum(sum(len_route, 2) == 0)
    length(id_leaves)

    %% Distance matrix

    leaves_names = arrayfun(@(x) {['L', num2str(x)]}, 1:n_leaves);
    len_route_leaves = len_route(id_leaves,:);
    points_leaves = []
    for i = 1:length(id_leaves)
        points_leaves(i,:) = segments_new(points(:,2) == id_leaves(i), 3:4);
    end

    writematrix(points_leaves, [path_tmp 'points_leaves.txt'])
    writematrix(len_route_leaves, [path_tmp 'len_route_leaves.txt'])


    % Plot tree and leaves + root
    f1 = figure; hold on;
    contour(X,Y,Z, 50);
    colormap(bone)
    for j = 1:size(segments_new, 1)
        plot([segments_new(j,1) segments_new(j,3)], [segments_new(j,2)  segments_new(j,4)], 'k-',...
            'LineWidth', 1.5)
    end
    axis equal
    plot(points_leaves(:,1), points_leaves(:,2), 'ok', 'MarkerFaceColor', 'k')
    plot(segments_new(1,1), segments_new(1,2), 'or', 'MarkerFaceColor', 'r')
%     text(points_leaves(:,1), points_leaves(:,2), leaves_names)
    axis equal
    plot(m(:,1), m(:,2), 'k+')
    saveas(gcf, [path_tmp 'segments_5_opt_leaves.pdf'])
    close(f1)

    % Distances on routes
    d_routes = zeros(n_leaves);
    for i = 1:n_leaves
        for j = 1:n_leaves
            if j <= i
                continue
            end
            idx_diff = (len_route_leaves(i,:) ~= len_route_leaves(j,:));
            d_tmp = sum(len_route_leaves(i,idx_diff)) + sum(len_route_leaves(j,idx_diff));
            d_routes(i,j) = d_tmp;
            d_routes(j,i) = d_tmp;
        end
    end
    d_routes_center = sum(len_route_leaves, 2)
    writematrix(d_routes_center, [path_tmp 'd_routes_center.txt'])
    writematrix(d_routes, [path_tmp 'd_routes.txt'])


    % Distances on lines
    d_lin_tmp = squareform(pdist([segments_new(1,1:2);points_leaves]));
    d_lin_center = d_lin_tmp(2:(n_leaves+1),1);
    d_lin = d_lin_tmp(2:(n_leaves+1),2:(n_leaves+1))
    writematrix(d_lin_center, [path_tmp 'd_lin_center.txt'])
    writematrix(d_lin, [path_tmp 'd_lin.txt'])


    % t_lin = seqneighjoin(d_lin, 'equivar', leaves_names)
    % t_routes = seqneighjoin(d_routes, 'equivar', leaves_names)
    % 
    % 
    % plot(t_lin)
    % plot(t_routes)

    % Correspondance between linear distances and routes
    f1 = figure; plot(d_lin(:), d_routes(:), 'o')
    axis equal
    saveas(gcf, [path_tmp 'corresp_routes_lin.pdf'])
    close(f1)

    %% Cov matrix


    cv_routes = (repmat(d_routes_center, 1, n_leaves) + repmat(d_routes_center', n_leaves, 1) - d_routes) / 2;
    cv_lin = (repmat(d_lin_center, 1, n_leaves) + repmat(d_lin_center', n_leaves, 1) - d_lin) / 2;

    writematrix(cv_routes, [path_tmp 'cv_routes.txt'])
    writematrix(cv_lin, [path_tmp 'cv_lin.txt'])



end


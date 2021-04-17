%%

n = 20;
time_max = 25;
tmp = randperm(time_max)
time_split = sort([2, tmp(1:(n-2))])

classes = repmat(1, n, 1)
px = repmat(0, n, 1)
py = repmat(0, n, 1)
f = repmat(0, n, 1)

segments = [];
scale = 1;
scalef = 0.5;

% rng(123)
for t = 2:time_max
    t
    classes(:,t) = classes(:,t-1);
    if sum(t == time_split) > 0
        % split
        classes_u = unique(classes(:,t));
        classes_u = classes_u(randperm(length(classes_u)));
        for cl = 1:length(classes_u)
            if sum(classes(:,t) == classes_u(cl)) > 0
                idx = 1:n;
                idx = idx(classes(:,t) == classes_u(cl));
                idx_bool = rand(length(idx),1) < 0.5;
                while (sum(idx_bool) == 0) || (sum(~idx_bool) == 0)
                    idx_bool = rand(length(idx),1) < 0.5;
                end
                idx = idx(idx_bool);
                classes(idx, t) = max(classes_u) + 1;
                break
            end
        end
    end
    n_cl = max(classes(:,t));
    [classes_u, ia, ib] = unique(classes(:,t));
    x_cl = px(ia,t-1);
    y_cl = py(ia,t-1);
    if(t > 2)
        x_cl_0 = px(ia,t-2);
        y_cl_0 = py(ia,t-2);
    end
    
    for icl = 1:n_cl
        icl;
        res = 1;
%         while(res == 1)
            
            dx = randn(1) / scale;
            dy = randn(1) / scale;
            segment_tmp = [x_cl(icl) y_cl(icl) x_cl(icl)+dx y_cl(icl)+dy];
            d_prev = [  x_cl(icl) - x_cl_0(icl) y_cl(icl) - y_cl_0(icl)];


%             if size(segments, 1) == 0
%                 res = 0;
%             end
            
            if (d_prev(1) * dx + d_prev(2) * dy) / norm([dx, dy]) / norm(d_prev) <= 2/3
                res = 1;
            end
            
            for iprev = 1:size(segments, 1)
                a = segment_tmp(1:2);
                b = segment_tmp(3:4);
                c = segments(:,1:2);
                d = segments(:,3:4);

                cross_prod = @(a1, a2)a1(:,1).*a2(:,2) - a1(:,2).*a2(:,1);

                % see c and d aroud (ab)
                ab = b - a;
                ac = c - a;
                ad = d - a;

                % see c and d aroud (cd)
                cb = b - c;
                ca = a - c;
                cd = d - c;

                if cross_prod(ac, ab) .* cross_prod(ad, ab) >= 0
                    res = 0;
                elseif cross_prod(ca, cd) .* cross_prod(cb, cd) >= 0
                    res = 0;
                else
                    res = 1;
                    break
                end

                figure; hold on;
                title(res)
                plot([a(1) b(1)], [a(2) b(2)])
                plot([c(1) d(1)], [c(2) d(2)])
            end
%         end
        if res == 1
            dx = 0
            dy = 0
        end
        x_cl(icl) = x_cl(icl)+dx;
        y_cl(icl) = y_cl(icl)+dy;
    end
    
    
    d_cl = sqrt(x_cl.^2 + y_cl.^2);
    f_cl = f(ia,t-1) + randn(n_cl,1) .* sqrt(d_cl) / scalef;
    
    
    px(:,t) = x_cl(ib);
    py(:,t) = y_cl(ib);
    f(:,t) = f_cl(ib);
    

    
    segments_tmp = [px(ia,t-1), py(ia,t-1), x_cl, y_cl];
    segments = [segments; segments_tmp];
    
    
    % add diffusion
end

plot(px', py')







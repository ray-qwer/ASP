function IMF_set = hht(x, t, thr)
    k = 0;
    k_max = 100;
    IMF_set = zeros(0, length(t));
    y = x;
    while true
        pre = [0,y(2:end)- y(1:end-1)];
        pro = [y(1:end-1)-y(2:end),0];
        local_max = find((pre > 0) & (pro > 0));
        local_min = find((pre < 0) & (pro < 0));
        if (length(local_max)+ length(local_min)<=1 || length(local_min)== 1 || length(local_max)==1)
            IMF_set = cat(1, IMF_set, y);
            break;
        end
        max_envelop = spline(t(local_max), y(local_max), t);
        min_envelop = spline(t(local_min), y(local_min), t);
        mean_line = (max_envelop+ min_envelop) ./2;
        h = y - mean_line;
        max_envelop = spline(t(local_max), h(local_max), t);
        min_envelop = spline(t(local_min), h(local_min), t);
        % check IMF
        if all(h(local_max) > 0) && all(h(local_min) < 0) && all(abs(max_envelop+ min_envelop)./2 < thr)
            % is IMF
            k = 0;
            c = h;
            IMF_set = cat(1, IMF_set, c);
            x = x - c;
            y = x;
        else
            k = k +1;
            y = h;
            if k > k_max 
                disp("the function is too weird");
                break;
            end
        end
    end
end
function [l_map_MI] = PCR_build_relaxed_edge_map(n,p_value,norm_MI,alpha)
    ss = xor(true(n), diag(true(1, n)));
    l_map_MI = zeros(0,3);
    l_cnt = 0;
    for i = 1:n-1
        for j = i+1:n
            if p_value(i, j) > alpha
                ss(i, j) = false; ss(j, i) = false;
            end
        end
    end
    for j = 1 : n-1
        for k = j+1 :n
            if ss(j,k)
                l_cnt = l_cnt + 1;
                l_map_MI(l_cnt,:) = [j,k,norm_MI(j,k)];
            end
        end
    end
end

function [p,l_map_MI,l_cnt] = CS_CGA_initialize_population(ss,N,norm_MI)
bns = size(ss,1);
p = cell(1,N);
p(:) = {false(bns)};
l_cnt = 0;
l_map_MI = zeros(0,3);
temp_a = rand(5*N,1);
temp_idx = find(temp_a>-0.5 & temp_a<0.5);
p_a = temp_a(temp_idx(1:N))+0.5;
for j = 1 : bns-1
    for k = j+1 :bns
        if ss(j,k)
            l_cnt = l_cnt + 1;
            l_map_MI(l_cnt,:) = [j,k,norm_MI(j,k)];
        end
    end
end
for l = 1:l_cnt
    j = l_map_MI(l,1);  k = l_map_MI(l,2);  l_MI = l_map_MI(l,3);
    for i = 1:N
        if l_MI > p_a(i)
            switch(randi(2))
                case 1
                    p{i}(j,k) = true;  p{i}(k,j) = false;
                case 2
                    p{i}(j,k) = false;  p{i}(k,j) = true;
                otherwise
                    p{i}(j,k) = false;  p{i}(k,j) = false;
            end
        end
    end
end
end
